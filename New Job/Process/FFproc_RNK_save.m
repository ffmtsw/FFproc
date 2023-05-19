function FFproc_RNK_save(final,stat)
    
    % Extracting information from final structure
    job = final.job;
    % Target frequencies
    target_freq = final.target_freq;
    time = final.time;

    % Number of time series per job
    nts = numel(job.ffts);

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

    % Number of oscillations
    osc = job.osc;

    RNK = final.RNK;
    
    % Magnetic reference channels
    magref = zeros(numel(job.ffts),2);
    for j = 1:numel(job.ffts)
        magref(j,1) = job.ffts(j).rsbx;
        magref(j,2) = job.ffts(j).rsby;
    end

    %% Saving estimations in corresponding files
    % Path for each FFproc user (Documents)
    load('ffprocsavepath.mat','path')
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
  
    % Saving Fourier Coefficients
    savename = [model,reg,num2str(osc),'OSC_',...
                num2str(f1),unit1,'-',num2str(f2),unit2,...
                remote,'_'...
                approach,'_',zrank,'_',trank,'_',evpow,'_',pcpow,...
                '_RNK100%'];

    % Erasing data for saving space
    for i = 1:nts
        job.ffts(i).t = [];
        job.ffts(i).data = [];
    end

    name = cat(1,job.name);
    save(savename,'RNK','job','target_freq','name','time','stat','-v7.3');

end