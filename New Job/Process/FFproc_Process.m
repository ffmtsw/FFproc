% function FFproc_Process executes the multi-variate remote reference
% processing for the submitted jobs created with FFproc_NewJob.
%
% Workflow:
%
%      - Loading time series
%      - Time series synchronization
%
%                      TIME DOMAIN
%           - Segmentation (time domain)
%               - Detrending (time domain)
%               - Filtering: Low/high pass (time domain)
%               - Tapering (time domain)
%
%                    FREQUENCY DOMAIN
%           - Fourier Transform (frequency domain)
%           - Calibration (frequency domain)
%           - Filtering: Unwanted frequencies (frequency domain)
%
%                    KERNEL FUNCTIONS
%           - Kernel 1: Multivariate Linear Regression (frequency domain)
%           - Kernel 2: Transfer Function Estimation (frequency domain)
%           - Kernel 3: Eigenvalue and Partial Coherence Evaluation (frequency domain)
%
%      - Merging results
%      - Final Transfer Functions: polyfit and sld
%
% FFproc is parallelized over time segments. Parallel pool is enabled by 
% default and used in Kernel 1, Kernel 2 and Kernel 3.
% 
% Indices used:
%      - j for stations
%      - k for channels (e.g. Ex, Ey, Bx, By, Bz)
%      - f for target frequencies
%      - w for time windows
%      - a for stations counter
%      - b for channel counter
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
%                                    FFproc
%                 An improved multivariate robust statistical 
%                 data processing software for the estimation 
%                          of MT transfer functions
%
%                                      .*(%,                               
%                      .%%%%%%%%%%%%%%%%%%%%%%.%%(                         
%                     %%%%%     %%%%%%%%%%%%%%%,(%%%%%%(                   
%                  .%%%%%%%%./#,  #%%. .%%%%  %%%%%%%%%%%                  
%                %%%%%%%   %%%%%%%%%%%%%%%%%*.*   .%%%%%%%                 
%               %%%%%, %%%%/                 %.,%%%%%%%%%#                 
%             ,%%%%%%%%%%%                    %%%%,  %%%%%%%               
%             *%%%%%%%%%%%                   /%%%%%%  %%%%%%,              
%              %%%%%%%%%%.                     %%%%%%%%%%%%%%              
%             %%%%%%%%%%.    .,.           ,,  *%%%%   %%%%%/              
%            .%%%%%%%%%%  %    , %%%   %%%(    %%%%  (%%%%%%               
%             ,%%%%%%%%%  ##( %/% %%  ,%%. %*% %#%%% .%%%  %%%             
%             %%.%%%%%%%.          %   %%,  ,%%* %%%%%%, .#%%%*            
%            /%%   (%%%            (   %%%      *%%%%%%%%%%%%%             
%             %%%%%%%%%  ,         .,  #/(%     /%%%%%(/%%%%.              
%                     ,%%%           %%%%.     .%%%%%%#                    
%                  /%%%%%%/     (,/%/ #%*,#   %%%%%%%%%%%%%                
%                %%%%%%%%%%         .   .%%%%%%%%%%%%%%%%%%%               
%                 %%%%%%%%   %*         .%%%%%%%% #%%%%%%%%                
%                  #%%%%%%    ,%%#     /%%%%%%%    %%%%%%#                 
%                    %%%%%%      .#%%%%%%%%.      /%%%%%(                  
%                (%%%% %%%%.#         %(         %%%%%%%%%%/               
%          *%%%%%%%%%%%,.%%# .%.    %%( %     /%%. %%%%%%%%%%%%%#          
%     .%%%%%%%%%%%%%%%%%% #%%    #%    (#%%%%%%  ,%%%%%%%%%%%%%%%%%%%#     
%  ,%%%%%%%%%%%%%%%%%%%%%%( %%%/  .% ,/.  %(    %%%%%%%%%%%%%%%%%%%%%%%%%  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               Hering, P., Castro, C., Hogg, C., Junge, A.
%                  Goethe-Universität Frankfurt am Main
%                     Institut für Geowissenschaften
%                               2020 - 2023
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% version 1.0 / 01feb2021 / cc          First estable release. FFproc works
%                                       with Metronix (ADU-07), Phoenix
%                                       (MTU-5A, V8), LEMI, Geolore and
%                                       Artificial Synthetic data produced
%                                       with TSGen - FFMT.
% version 2.0 / 05apr2021 / cc          Several modifications on time
%                                       series segmentation. The time
%                                       series are not segmented based on
%                                       the device anymore, but with an 
%                                       adaptive ts segmentation approach.
%                                       Phoenix MTU5-C devices can be 
%                                       processed. Different selected 
%                                       segments of time can be processed.
%
%%                              FFproc

function status = FFproc_Process(app,job,jobid,parproc,log)
    
    tic
    if app.procprogress.CancelRequested
        status = FFproc_cancelrequest(app);
        if ~isnan(status)
            return    
        end
    end
    
    if ~isempty(log)
        log.ProgresslogTextArea.Value(end+1,1) = {['Start date/time: ',char(datetime('now'))]};
        log.ProgresslogTextArea.Value(end+1,1) = {''};
        log.ProgresslogTextArea.Value(end+1,1) = {'--------------------------------------------------------------'};
        log.ProgresslogTextArea.Value(end+1,1) = {['                      Processing job ',num2str(jobid)]};
        log.ProgresslogTextArea.Value(end+1,1) = {'--------------------------------------------------------------'};
        log.ProgresslogTextArea.Value(end+1,1) = {''}; pause(0.1)
    end
    
    try
        %% -------------------------------------------- Loading Time Series
        app.procprogress.Message = 'Loading time series...';

        % Number of stations/time series
        nts = numel(job.ffts);

        % Preallocating data (structures)        
        data = struct;
        act_chan = cell(nts,1);

        % Loading time series from job and saving in data variable
        for j = 1:nts
            path = job.ffts(j).path;
            file = job.ffts(j).file;
            log.ProgresslogTextArea.Value(end+1,1) = {['Loading time series for site: ',job.ffts(j).name,' ...']}; pause(0.1)
           
            try
                app.procprogress.Message = ['Loading time series: ',job.ffts(j).file];
            catch
                app.procprogress.Message = ['Loading time series: ';job.ffts(j).file'];
            end
            switch job.ffts(j).device
                case 'Phoenix MTU5-C'
                    data.ts(j) = read_MTU5C(path,job.sr,true,job.ffts(j).tstart,job.ffts(j).tend);
                otherwise
                switch job.ffts(j).mode
                    case 'single'
                        data.ts(j) = FFproc_tsinput(path,file,job.ffts(j));
                    case 'merge'
                        tsmerge = [];
                        filemerge = strsplit(job.ffts(j).file,',');
                        pathmerge = strsplit(job.ffts(j).path,',');
                        for jj = 1:numel(filemerge)
                            tsmerge = [tsmerge,FFproc_tsinput(pathmerge{jj},filemerge{jj},job.ffts(j))]; %#ok<AGROW> 
                        end
                        data.ts(j) = ts_merge([],tsmerge,'FFproc');
                end
            end         
            log.ProgresslogTextArea.Value(end,1) = {['Loading time series for site: ',job.ffts(j).name,' - COMPLETE -']};  pause(0.1)             
            
            % Active channels configured in Site Information table
            act_chan{j} = job.ffts(j).actchan;
        end          
         
        if app.procprogress.CancelRequested
            status = FFproc_cancelrequest(app);
            if ~isnan(status)
                return    
            end
        end

        % Removing extra data in time series (make even number length)
        for i = 1:nts
            sizemax = floor(size(data.ts(i).data,1)/2)*2;
            data.ts(i).t = data.ts(i).t(1:sizemax);
            data.ts(i).data = data.ts(i).data(1:sizemax,:);
            data.ts(i).sn = sizemax;
        end      

        %% ------------------------------------------- Data synchronization
        app.procprogress.Message = 'Synchronizing time series...';  
        log.ProgresslogTextArea.Value(end+1,1) = {'Synchronization: ...'}; pause(0.1)

        if nts > 1
            % Number of workers
            if parproc         
                ncores = job.ncores;
            else
                ncores = 0;
            end
            
            [data.ts,msg] = ts_sync(data.ts,ncores);
            if msg == 2
                log.ProgresslogTextArea.Value(end+1,1) = {'Synchronization: Failed to syncrhonize!'}; pause(0.1)
                uialert(app.UIFigure,'At least one TS does not have common measured time',...
                                     ['FFproc - Synchronization Error in job: ',num2str(jobid)],'Icon','Error')
                status = 2;
                return
            elseif msg == 3
                log.ProgresslogTextArea.Value(end+1,1) = {'Synchronization: Failed to syncrhonize!'}; pause(0.1)
                uialert(app.UIFigure,'At least one TS presents problems - check for non-completed bursts',...
                                     ['FFproc - Synchronization Error in job: ',num2str(jobid)],'Icon','Error')
                status = 2;
                return
            else
                log.ProgresslogTextArea.Value(end,1) = {'Synchronization: - COMPLETE -'}; pause(0.1)
            end  
        else
            % Removing extra data out of specified start and end date and time
            for i = 1:nts
                % tstart
                tstartind = data.ts(i).t < job.ffts(i).tstart;
                data.ts(i).t(tstartind) = [];
                data.ts(i).data(tstartind,:) = [];
                data.ts(i).sn = numel(data.ts(i).t);
                clear tstartind
                
                % tend
                tendtind = data.ts(i).t > job.ffts(i).tend;
                data.ts(i).t(tendtind) = [];
                data.ts(i).data(tendtind,:) = [];
                data.ts(i).sn = numel(data.ts(i).t);
                clear tendtind
            end
            % Find gaps
            data.ts = findgaps(data.ts);            
        end
        if app.procprogress.CancelRequested
            status = FFproc_cancelrequest(app);
            if ~isnan(status)
                return    
            end
        end

        %% ---------------------------- Check for non-symmetric time series
        % It has been found that after the synchronization, some time
        % series has non-completed burst. This step tries to delete the
        % first and last burst
        data = ts_nonsymmetric(data,nts);
        
        %% ------------------------------------------- Time series schedule
        % Accounting for scheduled time series
        data = ts_schedule(data,nts,job);
                
        %% --------------------------------------------- Scaling E channels 
        % If the dipole lengths are modified from the original values,
        % there is a recalculation of the Electric field
        data = ts_echanscale(data,nts,job);
        
        %% ---- Accounting for errors on site deployment (swapped channels)
        data = ts_fixerrors(data,nts,job);

        %% Core Calculations following Egbert (1997) modified by Hering (2019)
        app.procprogress.Message = 'Kernel Calculations...';       
        
        % Definition of target frequencies
        [tfreq_lim,tfreq_vect,tfreq_cutoff] = find_tfreq(job.ftarg,job.maxf,job.minf);
        
        % Preallocating variables (vector, matrices and cell arrays)
        [n_tseg,n_samples,t_seg,osc,n_skip,overlap_tf,linreg,time,FC,TF,RNK,EigVal,results] = FFproc_allocatevars(tfreq_vect);
                       
        log.ProgresslogTextArea.Value(end+1,1) = {''};
        log.ProgresslogTextArea.Value(end+1,1) = {' ********** Start of transfer functions estimations ********* '}; pause(0.1)
        for f = 1:numel(tfreq_vect)
            % -------------------------------------- Target frequency label
            if tfreq_vect(f) > 1
                fqdisp = num2str(round(tfreq_vect(f),2));
                fqlabel = ' Hz';
            else
                fqdisp = num2str(round(1/tfreq_vect(f),2));
                fqlabel = ' s';
            end            
            log.ProgresslogTextArea.Value(end+1,1) = {['Estimations for target frequency: ',num2str(f),' of ',num2str(numel(tfreq_vect)),' @',fqdisp,fqlabel,' ...']}; pause(0.1)

            % ------------- Time series segmentation, detrending, filtering
            %               and tapering
            app.procprogress.Message = {'Segmentation, detrending, filtering and tapering';['for target frequency: ',num2str(f),' of ',num2str(numel(tfreq_vect)),' @',fqdisp,fqlabel]};
            if app.procprogress.CancelRequested
                status = FFproc_cancelrequest(app);
                if ~isnan(status)
                    return    
                end
            end
            [D,n_tseg(f),n_samples(f),t_seg(f),osc(f),n_skip(f),~,overlap_tf(f),time{f}] = ts_segment(data.ts,job.ffts,job.sr,tfreq_lim(f,:),tfreq_vect(f),tfreq_cutoff(f,:),f,job.osc,true); 
           
            % -------- Calculating spectral lines for each target frequency
            % Single-sided frequency vector calculated as fN*(1:N)/N
            freq = job.sr*(1:n_samples(f))/n_samples(f);   

            % Samples id of frequency vector used to estimate each target
            % Frequency between upper freq limit and Low-pass cutoff freq
            [~,freqmax] = min(abs(freq - 10^mean(log10([tfreq_cutoff(f,1),tfreq_lim(f,1)]))));            
            % Frequency between lower freq limit and Hi-pass cutoff freq
            [~,freqmin] = min(abs(freq - 10^mean(log10([tfreq_cutoff(f,2),tfreq_lim(f,2)]))));
            
            % ------------------ Segment of the frequency vector to be used
            freq_vec = freq(freqmin:freqmax); 

            % ------------- Fourier transform of each channel for each
            %               station
            clear F
            F = cell(nts,max(cat(1,data.ts(1).nchan)));
            for j = 1:nts
                app.procprogress.Message = {'Fourier Transform and Calibration';['for target frequency: ',num2str(f),' of ',num2str(numel(tfreq_vect)),' @',fqdisp,fqlabel]};
                if app.procprogress.CancelRequested
                    status = FFproc_cancelrequest(app);
                    if ~isnan(status)
                        return    
                    end
                end 

                % ------------------------ Tukey window for spectral domain
                fwin = repmat(tukeywin(numel(freq_vec),1),1,n_tseg(f));
                
                % --------------------------------------------- Calibration
                try
                    C = calibration(n_samples(f),data.ts(j).cal);                                       
                catch
                    C = calibration(n_samples(f),job.ffts(j).cal); 
                end

                % Calculation of the Fourier Transform for each time
                % segment (w), for each channel (k) and for each site (j)
                for k = 1:data.ts(j).nchan
                    % Extracting each channel with n_tseg time segments
                    dmat = D{j,k};     
                    % Fourier Transform of each data channel Dmat   
                    Dmat = fft(dmat);                    
                    % Calibration matrix
                    Cmat = repmat(C{k},1,n_tseg(f));                
                    % Matrix with calibrated Fourier Coefficients
                    Fmat = Dmat.*Cmat;
                    % Taking only 1/2 of the fft
                    Fmat = Fmat(1:end/2,:);
                    % Cell array with j sites and k channels. Each entry
                    % has only Fourier Coefficients around the target
                    % frequency.
                    F{j,k} = Fmat(freqmin:freqmax,:).*fwin;
                    % Releasing memory
                    clear dmat Fmat Dmat Cmat
                end
                
            end 
            clear C D
            % From now on, the units are:
            % - mV/km   for Electric field (E)
            % - nT      for Magnetic field (B)
            % - km/s    for Impedances (Z)
            % - Ohm*m   for App. Resistivities (rhoa)

            % ----------------------------------- Excluding Frequencies
            app.procprogress.Message = {'Excluding frequencies for target frequency: ';[num2str(f),' of ',num2str(numel(tfreq_vect)),' @',fqdisp,fqlabel]};
            if app.procprogress.CancelRequested
                status = FFproc_cancelrequest(app);
                if ~isnan(status)
                    return    
                end
            end    

            if isfield(job,'exclude')
                if job.exclude || isempty(job.exclude)        
                    remfq = unwantedfq(job,freq_vec,n_tseg(f));
                    for j = 1:nts
                        for k = 1:data.ts(j).nchan
                            F{j,k} = F{j,k}.*remfq;
                        end
                    end
                end
            end
            
            % -------------------- Kernel 1 Calculations: SDM, NCM, and Z
            app.procprogress.Message = {'Kernel 1: Multivariate Linear Regression (MVLR)';...
                                       ['calculation for target frequency: ',num2str(f),' of ',num2str(numel(tfreq_vect)),' @',fqdisp,fqlabel]};
            if app.procprogress.CancelRequested
                status = FFproc_cancelrequest(app);
                if ~isnan(status)
                    return    
                end
            end      
            [linreg{f},FC{f}] = FFproc_K1_ANM(job,F.',act_chan,n_tseg(f),parproc);
            clear F
            
            % -------------------- Kernel 2 Calculations: Transfer Function
            %                                             Estimations
            app.procprogress.Message = {'Kernel 2: Transfer Function Estimation';
                                       ['calculation for target frequency: ',num2str(f),' of ',num2str(numel(tfreq_vect)),' @',fqdisp,fqlabel]};
            if app.procprogress.CancelRequested
                status = FFproc_cancelrequest(app);
                if ~isnan(status)
                    return    
                end
            end
            [TF{f},FC{f}] = FFproc_K2_TF(job,linreg{f},act_chan,n_tseg(f),FC{f},parproc);

            % ---------------- --- Kernel 3 Calculations: Eigenvalue Criteria
            app.procprogress.Message = {'Kernel 3: Ranking Evaluation (EVI/PCOH)';...
                                       ['calculation for target frequency: ',num2str(f),' of ',num2str(numel(tfreq_vect)),' @',fqdisp,fqlabel]};
            if app.procprogress.CancelRequested
                status = FFproc_cancelrequest(app);
                if ~isnan(status)
                    return    
                end
            end
            [results{f},EigVal{f},RNK{f}] = FFproc_K3_RANK(job,n_tseg(f),act_chan,TF{f},RNK{f},parproc);

            % ----------------------- Calculations finished for Time Window
            log.ProgresslogTextArea.Value(end,1) = {['Estimations for target frequency: ',num2str(f),' of ',num2str(numel(tfreq_vect)),' @',fqdisp,fqlabel,' - COMPLETE -']}; pause(0.1)
            app.procprogress.Message = ['Calculation for window: ',num2str(f),' of ',num2str(numel(tfreq_vect)),' @',fqdisp,fqlabel,' has finished'];
            if app.procprogress.CancelRequested
                status = FFproc_cancelrequest(app);
                if ~isnan(status)
                    return    
                end
            end
        end
        
        %% ------------------------------------------ Finishing estimations
        log.ProgresslogTextArea.Value(end+1,1) = {' *********** End of transfer functions estimations ********** '}; pause(0.1)        
        app.procprogress.Title = ['FFproc: finishing job ',num2str(jobid)];
        
        %% ------------------------------------------------ Merging Results
        app.procprogress.Message = 'Merging results...';
        if app.procprogress.CancelRequested
            status = FFproc_cancelrequest(app);
            if ~isnan(status)
                return    
            end
        end
        % Creating Results variable (for transfer functions)
        log.ProgresslogTextArea.Value(end+1,1) = {''};
        log.ProgresslogTextArea.Value(end+1,1) = {'Merging estimated transfer functions: ...'}; pause(0.1)
        RAW = FFproc_merge_results(results);
        
        % Creating Results variable (for eigenvalues)
        EV = FFproc_merge_EV(TF,EigVal);
        log.ProgresslogTextArea.Value(end,1) = {'Merging estimated transfer functions: - COMPLETE -'}; pause(0.1)

        %% ------------------------------------------------ Saving FIT File      
        if job.polyfit
            try
                app.procprogress.Message = 'Fitting final transfer functions...';
                if app.procprogress.CancelRequested
                    status = FFproc_cancelrequest(app);
                    if ~isnan(status)
                        return    
                    end
                end
                log.ProgresslogTextArea.Value(end+1,1) = {'Fitting final transfer functions: ...'}; pause(0.1)
                    src = 'FIT';
                    final = FFproc_TFpolyfit(job,tfreq_vect,n_tseg,t_seg,act_chan,time,RNK,RAW,EV,TF,FC);
                log.ProgresslogTextArea.Value(end,1) = {'Fitting final transfer functions: - COMPLETE -'}; pause(0.1)
    
                app.procprogress.Message = 'Saving FIT file...';
                if app.procprogress.CancelRequested
                    status = FFproc_cancelrequest(app);
                    if ~isnan(status)
                        return    
                    end
                end
                log.ProgresslogTextArea.Value(end+1,1) = {'Saving FIT file: ...'}; pause(0.1)
                    FFproc_OPT_save(final,src);        
                log.ProgresslogTextArea.Value(end,1) = {'Saving FIT file: - COMPLETE -'}; pause(0.1) 
            catch
                log.ProgresslogTextArea.Value(end,1) = {'Saving FIT file: - FAILED -'}; pause(0.1) 
            end
        end
        
        %% ------------------------------------------------ Saving SLD File  
        app.procprogress.Message = 'Extracting selected transfer functions...';
        if app.procprogress.CancelRequested
            status = FFproc_cancelrequest(app);
            if ~isnan(status)
                return    
            end
        end  
        try
            log.ProgresslogTextArea.Value(end+1,1) = {'Extracting selected transfer functions: ...'}; pause(0.1)
                src = 'SLD';
                final = FFproc_TFwmedian(job,tfreq_vect,n_tseg,t_seg,act_chan,time,src,RNK,RAW,EV,TF,FC);
            log.ProgresslogTextArea.Value(end,1) = {'Extracting selected transfer functions: - COMPLETE -'}; pause(0.1)  
    
            app.procprogress.Message = 'Saving SLD file...';
            if app.procprogress.CancelRequested
                status = FFproc_cancelrequest(app);
                if ~isnan(status)
                    return    
                end
            end      
            log.ProgresslogTextArea.Value(end+1,1) = {'Saving SLD file: ...'}; pause(0.1)
                FFproc_OPT_save(final,src);        
            log.ProgresslogTextArea.Value(end,1) = {'Saving SLD file: - COMPLETE -'}; pause(0.1) 
        catch
            log.ProgresslogTextArea.Value(end,1) = {'Saving SLD file: - FAILED -'}; pause(0.1) 
        end
        
        %% ------------------------------------------------ Saving ALL File  
        app.procprogress.Message = 'Extracting final transfer functions...';
        if app.procprogress.CancelRequested
            status = FFproc_cancelrequest(app);
            if ~isnan(status)
                return    
            end
        end   
        try
            log.ProgresslogTextArea.Value(end+1,1) = {'Extracting final transfer functions: ...'}; pause(0.1)
                src = 'ALL';
                final = FFproc_TFwmedian(job,tfreq_vect,n_tseg,t_seg,act_chan,time,src,RNK,RAW,EV,TF,FC);
            log.ProgresslogTextArea.Value(end,1) = {'Extracting final transfer functions: - COMPLETE -'}; pause(0.1)
            
            app.procprogress.Message = 'Saving ALL file...';
            if app.procprogress.CancelRequested
                status = FFproc_cancelrequest(app);
                if ~isnan(status)
                    return    
                end
            end  
            log.ProgresslogTextArea.Value(end+1,1) = {'Saving ALL file: ...'}; pause(0.1)
                FFproc_OPT_save(final,src);        
            log.ProgresslogTextArea.Value(end,1) = {'Saving ALL file: - COMPLETE -'}; pause(0.1)       
        catch
            log.ProgresslogTextArea.Value(end,1) = {'Saving ALL file: - FAILED -'}; pause(0.1)  
        end

        %% ------------------------------------------------ Saving RNK File        
        if app.procprogress.CancelRequested
            status = FFproc_cancelrequest(app);
            if ~isnan(status)
                return    
            end
        end

        if job.rnk
            stat = FFproc_savestats(n_tseg,n_samples,t_seg,osc,n_skip,overlap_tf,time);
            try
                app.procprogress.Message = 'Saving RNK file...';
                log.ProgresslogTextArea.Value(end+1,1) = {'Saving RNK file: ...'}; pause(0.1)
                FFproc_RNK_save(final,stat);        
                log.ProgresslogTextArea.Value(end,1) = {'Saving RNK file: - COMPLETE -'}; pause(0.1)   
            catch
                log.ProgresslogTextArea.Value(end,1) = {'Saving RNK file: - FAILED -'}; pause(0.1)  
            end
        end
        
        %% --------------------------------------------- Saving RAW results        
        if app.procprogress.CancelRequested
            status = FFproc_cancelrequest(app);
            if ~isnan(status)
                return    
            end
        end

        if job.raw
            try
                app.procprogress.Message = 'Saving RAW file...';
                log.ProgresslogTextArea.Value(end+1,1) = {'Saving RAW file: ...'}; pause(0.1)
                FFproc_RAW_save(final);
                log.ProgresslogTextArea.Value(end,1) = {'Saving RAW file: - COMPLETE -'}; pause(0.1)
            catch
                log.ProgresslogTextArea.Value(end,1) = {'Saving RAW file: - FAILED -'}; pause(0.1)
            end
        end

        %% ------------------------------------------------ Saving FC File        
        if app.procprogress.CancelRequested
            status = FFproc_cancelrequest(app);
            if ~isnan(status)
                return    
            end
        end

        if job.fc
            try
                app.procprogress.Message = 'Saving FC file...';
                log.ProgresslogTextArea.Value(end+1,1) = {'Saving FC file: ...'}; pause(0.1)
                FFproc_FC_save(final);        
                log.ProgresslogTextArea.Value(end,1) = {'Saving FC file: - COMPLETE -'}; pause(0.1)   
            catch
                log.ProgresslogTextArea.Value(end,1) = {'Saving FC file: -FAILED -'}; pause(0.1)   
            end
        end

        %% ---------------------------------------------- Saving Parameters
        app.procprogress.Message = 'Saving PRM file...';
        if app.procprogress.CancelRequested
            status = FFproc_cancelrequest(app);
            if ~isnan(status)
                return    
            end
        end
        try
            log.ProgresslogTextArea.Value(end+1,1) = {'Saving PRM file...'}; pause(0.1)
            FFproc_PRM_save(final);
            log.ProgresslogTextArea.Value(end,1) = {'Saving PRM file: - COMPLETE -'}; pause(0.1)
        catch
            log.ProgresslogTextArea.Value(end,1) = {'Saving PRM file: - FAILED -'}; pause(0.1)
        end
        
        %% ------------------------------------------------- End of process
        finish = toc;
        log.ProgresslogTextArea.Value(end+1,1) = {''};
        log.ProgresslogTextArea.Value(end+1,1) = {['End date/time: ',char(datetime('now'))]};
        log.ProgresslogTextArea.Value(end+1,1) = {['Elapsed time: ',char(duration(0,0,finish))]};
        log.ProgresslogTextArea.Value(end+1,1) = {''};
        status = 1; 
    catch ME        
        report = getReport(ME);       
        uialert(app.UIFigure,report,['FFproc - Error in job: ',num2str(jobid)],'Icon','Error','Interpreter','html')
        status = -999;
        log.ProgresslogTextArea.Value(end+1,1) = {''};
        log.ProgresslogTextArea.Value(end+1,1) = {''};
        log.ProgresslogTextArea.Value(end+1,1) = {' ************************** ERROR! ************************* '};
        log.ProgresslogTextArea.Value(end+1,1) = {''};
        log.ProgresslogTextArea.Value(end+1,1) = {' Calculations can not be performed'};
        log.ProgresslogTextArea.Value(end+1,1) = {' Check for errors on the main window'};
        log.ProgresslogTextArea.Value(end+1,1) = {''};
        log.ProgresslogTextArea.Value(end+1,1) = {' ----------------------------------------------------------- '};
        log.ProgresslogTextArea.Value(end+1,1) = {['                         Job ',num2str(jobid),' Failed']};
        log.ProgresslogTextArea.Value(end+1,1) = {' ----------------------------------------------------------- '};
        log.ProgresslogTextArea.Value(end+1,1) = {''};
    end
    
end