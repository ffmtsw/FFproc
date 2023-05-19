function FFproc_OPT_save(final,src)
    
    % Extracting information from final structure
    job = final.job;
    % Target frequencies
    target_freq = final.target_freq;
    % 
    time = final.time;
    % Number of oscillations
    osc = job.osc;
    
    switch src
        case 'FIT'
            ranking = num2str(job.ranking);
            if numel(ranking) == 1
                ranking = ['00',ranking];
            elseif numel(ranking) == 2
                ranking = ['0',ranking];
            end
        case 'SLD'
            ranking = num2str(job.ranking);
            if numel(ranking) == 1
                ranking = ['00',ranking];
            elseif numel(ranking) == 2
                ranking = ['0',ranking];
            end
        case 'ALL'
            ranking = num2str(100);
    end

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

    FFsave = final.FFsave;
    EV = final.EV;
    id = final.id;
    cf = final.cf;
    
    % Number of time series per job
    nts = numel(job.ffts);
    
    % Magnetic reference channels
    magref = zeros(numel(job.ffts),2);
    for j = 1:numel(job.ffts)
        magref(j,1) = job.ffts(j).rsbx;
        magref(j,2) = job.ffts(j).rsby;
    end

    %% Saving estimations in corresponding files
    % Path for each FFproc user (Documents)
    try
        load('ffprocsavepath.mat','path')
    catch ME    
    end
    cd(path)
    
    % Creating folder name for each job an creating folder (if it does not
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
    
    % Retrieving target frequency range label: seconds or Hz
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
    
    % Assembling output FFsave file name
    savename = [model,reg,num2str(osc),'OSC_',...
                num2str(f1),unit1,'-',num2str(f2),unit2,...
                remote,'_',...
                approach,'_',zrank,'_',trank,'_',evpow,'_',pcpow,'_'...
                src,ranking,'%'];
    
    % Erasing data for saving space
    for i = 1:nts
        job.ffts(i).t = [];
        job.ffts(i).data = [];
    end
    
    % Saving estimations
    switch src
        case 'ALL'
            save(savename,'job','FFsave','EV','id','cf','time','-v7.3');
        otherwise
            save(savename,'job','FFsave','id','cf','time','-v7.3');
    end

end