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
%                              FFproc

function status = FFproc_Process(~,job,~,parproc,log)
    
    tic    
    try
        % -------------------------------------------- Loading Time Series
        % Number of stations/time series
        nts = numel(job.ffts);

        % Preallocating data (structures)        
        data = struct;
        act_chan = cell(nts,1);

        % Loading time series from job and saving in data variable
        for j = 1:nts
            path = job.ffts(j).path;
            file = job.ffts(j).file;

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
            
            % Active channels configured in Site Information table
            act_chan{j} = job.ffts(j).actchan;
        end          
         
        % Removing extra data in time series (make even number length)
        for i = 1:nts
            sizemax = floor(size(data.ts(i).data,1)/2)*2;
            data.ts(i).t = data.ts(i).t(1:sizemax);
            data.ts(i).data = data.ts(i).data(1:sizemax,:);
            data.ts(i).sn = sizemax;
        end      

        % ------------------------------------------- Data synchronization
        if nts > 1
            % Number of workers
            if parproc         
                ncores = job.ncores;
            else
                ncores = 0;
            end
            
            [data.ts,msg] = ts_sync(data.ts,ncores);
            if msg == 2
                return
            elseif msg == 3
                status = 2;
                return
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

        % ---------------------------- Check for non-symmetric time series
        % It has been found that after the synchronization, some time
        % series has non-completed burst. This step tries to delete the
        % first and last burst
        data = ts_nonsymmetric(data,nts);
        
        % ------------------------------------------- Time series schedule
        % Accounting for scheduled time series
        data = ts_schedule(data,nts,job);
                
        % --------------------------------------------- Scaling E channels 
        % If the dipole lengths are modified from the original values,
        % there is a recalculation of the Electric field
        data = ts_echanscale(data,nts,job);
        
        % ---- Accounting for errors on site deployment (swapped channels)
        data = ts_fixerrors(data,nts,job);
        
        % Definition of target frequencies
        tfreq_lim = find_tfreq(job.ftarg,job.minf,job.maxf);
        tfreq_vect = mean(tfreq_lim,2);
        
        % Preallocating variables (vector, matrices and cell arrays)
        [n_tseg,n_samples,t_seg,linreg,time,FC,TF,RNK,EigVal,results] = FFproc_allocatevars(tfreq_vect);
                       
        % *************************************************************** %
        %                       START OF ESTIMATIONS                      %
        % *************************************************************** %
        for f = 1:numel(tfreq_vect)

            % Time series segmentation, detrending, filtering and tapering
            [D,n_tseg(f),n_samples(f),t_seg(f),~,~,~,~,time{f},cutoff] = ts_segment(data.ts,job.sr,tfreq_lim(f,:),tfreq_vect(f),f,job.osc,true); 

            % Spectral lines for each target frequency
            % (single-sided frequency vector calculated as fN*(1:N)/N)
            freq = job.sr*(1:n_samples(f))/n_samples(f);   

            % Samples id of frequency vector used to estimate each target 
            % Frequency between lower freq limit and Hi-pass cutoff freq
            [~,freqmin] = min(abs(freq - 10^mean(log10([cutoff(1),tfreq_lim(f,2)])))); 
            % Frequency between upper freq limit and Low-pass cutoff freq
            [~,freqmax] = min(abs(freq - 10^mean(log10([cutoff(2),tfreq_lim(f,1)]))));
            
            % Segment of the frequency vector to be used
            freq_vec = freq(freqmin:freqmax); 

            % Fourier transform of each channel for each station
            clear F
            F = cell(nts,max(cat(1,data.ts(1).nchan)));
            for j = 1:nts
                
                % Calibration
                try
                    C = calibration(n_samples(f),data.ts(j).cal);                                       
                catch
                    C = calibration(n_samples(f),job.ffts(j).cal); 
                end

                % Calculate the Fourier Transform for each time segment 
                % (w), for each channel (k) and for each site (j)
                for k = 1:data.ts(j).nchan
                    % Extract each channel with n_tseg time segments
                    dmat = D{j,k};     
                    % Fourier Transform of each data channel Dmat   
                    Dmat = fft(dmat);                    
                    % Calibration matrix
                    Cmat = repmat(C{k},1,n_tseg(f));                
                    % Matrix with calibrated Fourier Coefficients
                    Fmat = Dmat.*Cmat;
                    % Take only 1/2 of the fft
                    Fmat = Fmat(1:end/2,:);
                    % Cell array with j sites and k channels. Each entry
                    % has only Fourier Coefficients around the target
                    % frequency.
                    F{j,k} = Fmat(freqmin:freqmax,:);
                    % Release memory
                    clear dmat Fmat Dmat Cmat
                end
                
            end 
            clear C D
            % From now on, the units are:
            % - mV/km   for Electric field (E)
            % - nT      for Magnetic field (B)
            % - km/s    for Impedances (Z)
            % - Ohm*m   for App. Resistivities (rhoa)

            % Exclude unwanted frequencies
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
            
            % Kernel 1 Calculations: SDM, NCM, and Z
            [linreg{f},FC{f}] = FFproc_K1_ANM(job,F.',act_chan,n_tseg(f),parproc);
            clear F            
            % Kernel 2 Calculations: Transfer Function Estimations
            [TF{f},FC{f}] = FFproc_K2_TF(job,linreg{f},act_chan,n_tseg(f),FC{f},parproc);
            % Kernel 3 Calculations: Eigenvalue Criteria
            [results{f},EigVal{f},RNK{f}] = FFproc_K3_RANK(job,n_tseg(f),act_chan,TF{f},RNK{f},parproc);
            
        end
        
        % *************************************************************** %
        %                        END OF ESTIMATIONS                       %
        % *************************************************************** %
        
        %  Merge Results
        % Creating Results variable (for transfer functions)
        RAW = FFproc_merge_results(results);        
        % Create Results variable (for eigenvalues)
        EV = FFproc_merge_EV(TF,EigVal);

        % ------------------------------------------------ Saving FIT File      
        src = 'FIT';
        if job.polyfit
            try                
                final = FFproc_TFpolyfit(job,tfreq_vect,n_tseg,t_seg,act_chan,time,RNK,RAW,EV,TF,FC);
                FFproc_OPT_save(final,src);        
            catch
            end
        end
        
        % ------------------------------------------------ Saving SLD File  
        src = 'SLD';
        try            
            final = FFproc_TFwmedian(job,tfreq_vect,n_tseg,t_seg,act_chan,time,src,RNK,RAW,EV,TF,FC); 
            FFproc_OPT_save(final,src);        
        catch
        end
        
        % --------------------------- ALL File -------------------------- % 
        src = 'ALL';
        try           
            final = FFproc_TFwmedian(job,tfreq_vect,n_tseg,t_seg,act_chan,time,src,RNK,RAW,EV,TF,FC);
            FFproc_OPT_save(final,src);        
        catch
        end

        % --------------------------- RNK File -------------------------- %       
        if job.rnk
            try FFproc_RNK_save(final); catch; end
        end
        
        % --------------------------- RAW File -------------------------- %       
        if job.raw
            try FFproc_RAW_save(final); catch; end
        end

        % --------------------------- FC File --------------------------- %
        if job.fc
            try FFproc_FC_save(final); catch; end
        end

        % --------------------------- PRM File -------------------------- %
        try FFproc_PRM_save(final); catch; end
        
        % ------------------------------------------------- End of process
        finish = toc;
        log.ProgresslogTextArea.Value(end+1,1) = {''};
        disp(['End date/time: ',char(datetime('now'))]);
        disp(['Elapsed time: ',char(duration(0,0,finish))]);
        status = 1; 

    catch ME        
        report = getReport(ME);       
        status = -999;
    end
    
end