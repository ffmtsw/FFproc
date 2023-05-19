function FFproc_merge_RAW(app,file,path)

    % Loading files and retrieving information 
    for i = 1:numel(path)
        load(fullfile(path{i},file{i}))
        job_merge{i} = job;
        act_chan_merge{i} = act_chan; 
        RAW_merge{i} = RAW;
        EV_merge{i} = EV;        
    end
    for i = 1:numel(EV_merge)
        sr(i) = EV_merge{i}(1).sr;
    end
    [~,sr_id] = sort(sr,'descend');
    
    % Arranging files
    RAW_merge = RAW_merge(sr_id);
    EV_merge = EV_merge(sr_id);
    
    % Ammount of best ranked EV
    for i = 1:numel(EV_merge)
        perc(i) = 100;
    end
    
    if numel(unique(perc)) ~= 1
        answer = uiconfirm(app.UIFigure,{'RAW Files contain different amount of ranked EV';...
                                         'Only common ranked EV will be taken into account';...
                                         'Do you want to proceed?'}, ...
                                         'FFproc - Merge RAW Files', ...
                                         'Options',{'Yes','No'},'DefaultOption',1,'CancelOption',2);
        switch answer
            case 'Yes'
                perc = min(perc);
            case 'No'
                return
        end
    end
    
    % Cleaning memory
    clear job act_chan RAW 
    
    % Concatenating RAW structure
    RAW = create_ffsave;
    for j = 1:perc
        
        % Allocating variables
        site = [];      lonlat = [];    UTM = [];       z = [];
        nfreq = [];     freq = [];      per = [];
        Zxx = [];       Zxy = [];       Zyx = [];       Zyy = [];
        Zxx_Err = [];   Zxy_Err = [];   Zyx_Err = [];   Zyy_Err = [];
        txz = [];       txz_Err = [];
        tyz = [];       tyz_Err = [];   
        txx = [];       txy = [];       tyx = [];       tyy = [];
        t_id = [];
        phi11 = [];     phi12 = [];     phi21 = [];     phi22 = [];        
        alpha = [];     beta = [];      theta = [];     lambda = [];
        alpha_Err = []; beta_Err = [];  theta_Err = []; lambda_Err = [];
        phimax = [];    phimin = [];        
        phimax_Err = [];phimin_Err = [];        
        
        for i = 1:numel(RAW_merge)
            
            % Name and coordinates
            site = [site;RAW_merge{i}(j).site];
            lonlat = [lonlat;{RAW_merge{i}(j).lonlat}];
            UTM = [UTM;{RAW_merge{i}(j).UTM}];
            z = [z;{RAW_merge{i}(j).z}];
            
            % Frequencies
            nfreq = [nfreq;RAW_merge{i}(j).nfreq];
            freq = [freq;RAW_merge{i}(j).freq];
            per = [per;RAW_merge{i}(j).per];
            
            % Impedances
            Zxx = [Zxx;RAW_merge{i}(j).Zxx];
            Zxy = [Zxy;RAW_merge{i}(j).Zxy];
            Zyx = [Zyx;RAW_merge{i}(j).Zyx];
            Zyy = [Zyy;RAW_merge{i}(j).Zyy];
            Zxx_Err = [Zxx_Err;RAW_merge{i}(j).Zxx_Err];
            Zxy_Err = [Zxy_Err;RAW_merge{i}(j).Zxy_Err];
            Zyx_Err = [Zyx_Err;RAW_merge{i}(j).Zyx_Err];
            Zyy_Err = [Zyy_Err;RAW_merge{i}(j).Zyy_Err];
            
            % Tipper
            txz = [txz;RAW_merge{i}(j).txz];
            tyz = [tyz;RAW_merge{i}(j).tyz];
            txz_Err = [txz_Err;RAW_merge{i}(j).txz_Err];
            tyz_Err = [tyz_Err;RAW_merge{i}(j).tyz_Err];
            
            % Phase Tensor
            phi11 = [phi11;RAW_merge{i}(j).phi11];
            phi12 = [phi12;RAW_merge{i}(j).phi12];
            phi21 = [phi21;RAW_merge{i}(j).phi21];
            phi22 = [phi22;RAW_merge{i}(j).phi22];
            alpha = [alpha;RAW_merge{i}(j).alpha];
            beta = [beta;RAW_merge{i}(j).beta];
            theta = [theta;RAW_merge{i}(j).theta];
            lambda = [lambda;RAW_merge{i}(j).lambda];
            phimax = [phimax;RAW_merge{i}(j).phimax];
            phimin = [phimin;RAW_merge{i}(j).phimin];
            alpha_Err = [alpha_Err;RAW_merge{i}(j).alpha_Err];
            beta_Err = [beta_Err;RAW_merge{i}(j).beta_Err];
            theta_Err = [theta_Err;RAW_merge{i}(j).theta_Err];
            lambda_Err = [lambda_Err;RAW_merge{i}(j).lambda_Err];
            phimax_Err = [phimax_Err;RAW_merge{i}(j).phimax_Err];
            phimin_Err = [phimin_Err;RAW_merge{i}(j).phimin_Err];
            
            % HMTF
            txx = [txx;RAW_merge{i}(j).txx];       
            txy = [txy;RAW_merge{i}(j).txy]; 
            tyx = [tyx;RAW_merge{i}(j).tyx];         
            tyy = [tyy;RAW_merge{i}(j).tyy]; 
            t_id = [t_id;RAW_merge{i}(j).t_id];
            
        end
        
        % Name and coordinates
        RAW(j).site = unique(site);
        RAW(j).nfreq = unique(nfreq);
        RAW(j).lonlat = lonlat{1};
        RAW(j).UTM = UTM{1};
        RAW(j).z = z{1};
        
        % Frequencies
        RAW(j).nfreq = sum(nfreq);
        RAW(j).freq = freq;
        RAW(j).per = per;
        
        % Saving transfer functions in RAW structure: Z
        RAW(j).Zxx = Zxx;                           RAW(j).Zxx_Err = Zxx_Err;                       
        RAW(j).Zxy = Zxy;                           RAW(j).Zxy_Err = Zxy_Err;
        RAW(j).Zyx = Zyx;                           RAW(j).Zyx_Err = Zyx_Err;        
        RAW(j).Zyy = Zyy;                           RAW(j).Zyy_Err = Zyy_Err;   
        
        % Saving transfer functions in RAW structure: T
        RAW(j).txz = txz;                            RAW(j).txz_Err = txz_Err;
        RAW(j).tyz = tyz;                            RAW(j).tyz_Err = tyz_Err;
        
        % Saving transfer functions in RAW structure: PT
        RAW(j).phi11 = phi11;                        RAW(j).phi12 = phi12;
        RAW(j).phi21 = phi21;                        RAW(j).phi22 = phi22;
        RAW(j).beta = beta;                          RAW(j).beta_Err = beta_Err;
        RAW(j).alpha = alpha;                        RAW(j).alpha_Err = alpha_Err;
        RAW(j).theta = theta;                        RAW(j).theta_Err = theta_Err;
        RAW(j).lambda = lambda;                      RAW(j).lambda_Err = lambda_Err;
        RAW(j).phimax = phimax;                      RAW(j).phimax_Err = phimax_Err;
        RAW(j).phimin = phimin;                      RAW(j).phimin_Err = phimin_Err;   
        
    end
    
    % Resistivity and phase calculations, scaling and corrections
    RAW = FFproc_RAW_rhoandphi(RAW,job_merge{1});   
    
%     clear EV
%     % Concatenating EV structure
%     for j = 2:perc
%         best = [];      
%         for i = 1:numel(RAW_merge)
%             % Best estimations
%             best = [best;EV_merge{i}.best(j).best];
%         end
%         EV.best(j).best = best;
%     end
%     
%     % Allocating variables
%     perc = [];      EVal = [];      Ev = []; 
%     timeseg = [];   tfreq = [];     ntseg = []; 
%     tstart = [];    tend = [];  
%     
%     for i = 1:numel(RAW_merge)
%         
%         perc = [perc,EV_merge{i}.perc];
%         EVal = [EVal,EV_merge{i}.EVal]; 
%         Ev = [Ev,EV_merge{i}.EV]; 
%         timeseg = [timeseg,EV_merge{i}.timeseg]; 
%         tfreq = [tfreq,EV_merge{i}.tfreq]; 
%         ntseg = [ntseg,EV_merge{i}.tfreq];
%         tstart = [tstart;EV_merge{i}.tstart];
%         tend = [tend;EV_merge{i}.tend];
%         
%     end
%     
%     EV.perc = min(perc);
%     EV.EVal = EVal;
%     EV.EV = EV;
%     EV.timeseg = timeseg;
%     EV.tfreq = tfreq;
%     EV.ntseg = ntseg;
%     EV.tstart = tstart;
%     EV.tend = tend;      
    
    % Saving magnetic reference status (to identify local or RR)
    for i = 1:numel(job_merge)
        for j = 1:numel(job_merge{i}.ffts)
            magref{i,j} = [job_merge{i}.ffts(j).rsbx,job_merge{i}.ffts(j).rsby];
            device{i,j} = job_merge{i}.ffts(j).device;
        end
    end  
    
    site = unique(site);    
    for i = 1:size(magref,2)
        mag(i,:) = all(cell2mat(magref(:,i)));     
        job.ffts(i).rsbx = mag(i,1);
        job.ffts(i).rsby = mag(i,2);
        job.name(i) = site(i);
        job.ffts(i).GPS = [lonlat(i,:),z(i)];
        job.ffts(i).device = char(unique(device));
        job.eigperc = min(perc);
        job.ftarg = 6;
        job.ffts(i).t = [];
        job.ffts(i).data = [];
    end
    
    % Saving active channels array (to identify which channels were used
    % to perform calculations
    for i = 1:numel(job_merge)
        for j = 1:numel(job_merge{i}.ffts)
            act{i,j} = act_chan_merge{1,i};
        end
    end
    for i = 1:size(act,2)
        act_chan{i} = all(reshape(cell2mat([act{:,i}]),5,[])');
    end    
    
    % Retrieving target frequency range label: seconds or Hz
    if max(freq) >= 1
        f1 = round(max(freq));
        unit1 = 'Hz';
    else
        f1 = round(1./max(freq));
        unit1 = 's';
    end
    if min(freq) >= 1
        f2 = round(min(freq));
        unit2 = 'Hz';
    else
        f2 = round(1./min(freq));
        unit2 = 's';
    end

    switch job_merge{1}.model
        case 'Standard'
            model = 'STD_';
        case 'Advanced'
            model = 'ANM_';
    end

    switch job_merge{1}.reg
        case 'Robust'
            reg = 'ROB_';
        case 'LSQ'
            reg = 'LSQ_';
    end
    
    % Saving file
    outname = ['MERGE_',model,reg,num2str(f1),unit1,'-',num2str(f2),unit2,'_RAW',num2str(min(perc)),'%'];       
    [file,path] = uiputfile([outname,'.mat'],'Save merged RAW *.mat file');
    cd(path)
    Results = RAW;
    target_freq = freq;
    EV = [];

    save(fullfile(path,file),'Results','EV','job','act_chan','target_freq');    
    
end