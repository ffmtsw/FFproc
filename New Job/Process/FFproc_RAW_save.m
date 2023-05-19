function FFproc_RAW_save(final)
    
    % Extracting information from final structure
    job = final.job;
    % Target frequencies
    target_freq = final.target_freq;
    % EigenValues
    EV = final.EV;
    % RAW estimations
    Results = final.RAW;
    % Active channels used for processing
    act_chan = final.act_chan;

    % Processing approach
    if job.bivar && ~job.ev         % Bivariate
        approach = 'BVR';
    elseif ~job.bivar && job.ev     % Eigenvalue Decomposition
        approach = 'EVD';
        approach = [approach,'-',job.eind];
    end

    % Transfer functions ranking
    zrank = ['ZRNK-',job.zrank];
    trank = ['TRNK-',job.trank];

    % Powers
    evpow = ['EIPOW-',num2str(job.evpow)];
    pcpow = ['PCPOW-',num2str(job.rpow)];

    % Number of time series per job
    nts = numel(job.ffts);    
    % Maximum of ammount of data to be considered for estimations (1-100%)
    ranking = num2str(numel(Results));
    if numel(ranking) == 1
        ranking = ['00',ranking];
    elseif numel(ranking) == 2
        ranking = ['0',ranking];
    end
    % Number of oscillations
    osc = job.osc;
    
    % Magnetic reference channels
    magref = zeros(numel(job.ffts),2);
    for j = 1:numel(job.ffts)
        magref(j,1) = job.ffts(j).rsbx;
        magref(j,2) = job.ffts(j).rsby;
    end
    
    % Create FFsave and EV structures to pre-allocate fields.
    RAW = create_ffsave;
    
    for r = 1:numel(Results)        
        % ------------------------------ Saving all the estimations for Zxx
        Zxx = Results(r).Zxx;
        Zxx_Err = Results(r).Zxx_Err;
        % ------------------------------ Saving all the estimations for Zxy
        Zxy = Results(r).Zxy;  
        Zxy_Err = Results(r).Zxy_Err;
        % ------------------------------ Saving all the estimations for Zyx
        Zyx = Results(r).Zyx;
        Zyx_Err = Results(r).Zyx_Err;
        % ------------------------------ Saving all the estimations for Zyy
        Zyy = Results(r).Zyy;     
        Zyy_Err = Results(r).Zyy_Err;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                PHASE TENSOR ESTIMATION                    %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        phi11 = Results(r).phi11;      phi12 = Results(r).phi12;
        phi21 = Results(r).phi21;      phi22 = Results(r).phi22;
        beta = Results(r).beta;        beta_Err = Results(r).beta_Err;
        alpha = Results(r).alpha;      alpha_Err = Results(r).alpha_Err;
        theta = Results(r).theta;      theta_Err = Results(r).theta_Err;
        lambda = Results(r).lambda;    lambda_Err = Results(r).lambda_Err;
        phimax = Results(r).phimax;    phimax_Err = Results(r).phimax_Err;
        phimin = Results(r).phimin;    phimin_Err = Results(r).phimin_Err;     

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                   TIPPER ESTIMATION                       %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Deleting previous estimations and preallocating variables
        clear txz tyz txz_Err tyz_Err        
        txz = Results(r).txz;                           txz_Err = Results(r).txz_Err;
        tyz = Results(r).tyz;                           tyz_Err = Results(r).tyz_Err;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                 SAVING FINAL ESTIMATIONS                  %%%%
        %%%%                 FOR EACH BEST ESTIMATION                  %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        % Saving number of frequencies, frequency and period vector
        RAW(r).nfreq = numel(target_freq);
        RAW(r).freq = target_freq;
        RAW(r).per = 1./(target_freq);
        % Saving transfer functions in FFsave structure: Z
        RAW(r).Zxx = Zxx;                            RAW(r).Zxx_Err = Zxx_Err;
        RAW(r).Zxy = Zxy;                            RAW(r).Zxy_Err = Zxy_Err;
        RAW(r).Zyx = Zyx;                            RAW(r).Zyx_Err = Zyx_Err;
        RAW(r).Zyy = Zyy;                            RAW(r).Zyy_Err = Zyy_Err;
        % Saving transfer functions in FFsave structure: T
        RAW(r).txz = txz;                            RAW(r).txz_Err = txz_Err;
        RAW(r).tyz = tyz;                            RAW(r).tyz_Err = tyz_Err;
        % Saving transfer functions in results structure: PT
        RAW(r).phi11 = phi11;                        RAW(r).phi12 = phi12;
        RAW(r).phi21 = phi21;                        RAW(r).phi22 = phi22;
        RAW(r).beta = beta;                          RAW(r).beta_Err = beta_Err;
        RAW(r).alpha = alpha;                        RAW(r).alpha_Err = alpha_Err;
        RAW(r).theta = theta;                        RAW(r).theta_Err = theta_Err;
        RAW(r).lambda = lambda;                      RAW(r).lambda_Err = lambda_Err;
        RAW(r).phimax = phimax;                      RAW(r).phimax_Err = phimax_Err;
        RAW(r).phimin = phimin;                      RAW(r).phimin_Err = phimin_Err;           
    end
    
    % Writing extra info to the resultant RAW variable:
    % Name, coordinates (lonlat, UTM) and elevation.
    RAW = FFproc_RAW_metainfo(RAW,job);
    
    % Resistivity and phase calculations, escaling and corrections
    RAW = FFproc_RAW_rhoandphi(RAW,job);
        
    %% Saving estimations in corresponding files
    % Path for each FFproc user (Documents)
    load('ffprocsavepath.mat','path')
    cd(path)
    
    % Creating folder name for each job an creating folder (if does not
    % exists)
    file = job.name;
    foldername = [];
    if numel(file) == 1
        foldername = file{1};
    else
        for n = 1:numel(job.name)
            foldername = [foldername,[file{n}]]; %#ok<*AGROW>
            if n < numel(file)
                foldername = [foldername,'-'];
            end
        end
    end
    if isfolder(foldername)
    else
        mkdir(foldername)
    end
    cd(foldername)
    
    % Retrieving target frequency range names
    if max(target_freq) >= 1
        f1 = round(target_freq(1));
        unit1 = 'Hz';
    else
        f1 = round(1./target_freq(1));
        unit1 = 's';
    end
    if min(target_freq) >= 1
        f2 = round(target_freq(end));
        unit2 = 'Hz';
    else
        f2 = round(1./target_freq(end));
        unit2 = 's';
    end
    
    % Remote Reference indicator
    if  numel(magref(:,1)) ~= numel(unique(magref(:,1))) || ...
        numel(magref(:,2)) ~= numel(unique(magref(:,2)))
        remote = '_RR';
    else
        remote = '';
    end
    
    switch job.model
        case 'Standard'
            model = 'STD_';
        case 'Advanced'
            model = 'ANM_';
    end

    switch job.reg
        case 'Robust'
            reg = 'ROB_';
        case 'LSQ'
            reg = 'LSQ_';
    end
    
    % Assembling output RAW file name
    savename = [model,reg,num2str(osc),'OSC_',...
                num2str(f1),unit1,'-',num2str(f2),unit2,...
                remote,'_'...
                approach,'_',zrank,'_',trank,'_',evpow,'_',pcpow,...
                '_RAW',ranking,'%'];
    
    % Saving estimations
    save(savename,'job','RAW','target_freq','act_chan','-v7.3');

end