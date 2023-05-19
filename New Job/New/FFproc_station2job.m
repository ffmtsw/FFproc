function FFproc_station2job(app,job,ffts,schedule)

    if isempty(job)
        job = job_create;
        nj = 1;
    else
        nj = numel(job) + 1;
    end

    % Station name
    name = cell(1,numel(ffts));
    for i = 1:numel(ffts)
        name{i} = ffts(i).name;
    end
    
    % From Target Frequencies
    maxf = app.maxfreqbox.Value;
    minf = app.minfreqbox.Value;
    ftarg = app.targfreqbox.Value;
    osc = app.oscbox.Value;

    % From Noise Removal
    exclude = app.exclude.State;
    switch exclude
        case 'on'
            excfreq = true;
        case 'off'
            excfreq = false;
    end
    mult = str2num(app.multbox.Value);
    singf = str2num(app.singfbox.Value);
    phmult = str2double(app.phmultbox.Value);
    phsing = str2double(app.phsingbox.Value);

    % From TF Estimation Approach
    if app.BivariateLinearRegressionButton.Value
        ev = false;
        bivar = true;
    elseif app.EigenvalueDecompositionButton.Value
        ev = true;
        bivar = false;
    end

    % From Model Estimator 
    if app.NMStandardButton.Value
        model = 'Standard';
    else
        model = 'Advanced';
    end      

    % From Linear Regression
    if app.LRLSQButton.Value
        reg = 'LSQ';
    elseif app.LRRobustButton.Value
        reg = 'Robust';
    elseif app.LRPCAButton.Value
        reg = 'PCA';
    end

    % From Outlier Detection
    outlier = logical(app.OutlierButton.Value);
    remmethod = app.remmethod.Value;

    % From Ranking Criterion
    if app.EVZButton.Value && ~app.PCCZButton.Value
        zrank = 'EV';
    elseif ~app.EVZButton.Value && app.PCCZButton.Value
        zrank = 'PCC';
    elseif app.EVZButton.Value && app.PCCZButton.Value
        zrank = 'EVPCC';
    end
    if app.EVTButton.Value && ~app.PCCTButton.Value
        trank = 'EV';
    elseif ~app.EVTButton.Value && app.PCCTButton.Value
        trank = 'PCC';
    elseif app.EVTButton.Value && app.PCCTButton.Value
        trank = 'EVPCC';
    end
    ranking = app.ranking.Value;
    
    % From Eigenvalue Index
    if app.EI0Button.Value
        eind0 = 'EI0';
    else
        eind0 = '';
    end
    if app.EI1Button.Value
        eind1 = 'EI1';
    else
        eind1 = '';
    end
    if app.EI2Button.Value
        eind2 = 'EI2';
    else
        eind2 = '';
    end
    if app.EI3Button.Value
        eind3 = 'EI3';
    else
        eind3 = '';
    end
    if app.EI4Button.Value
        eind4 = 'EI4';
    else
        eind4 = '';
    end

    eind = [eind0,eind1,eind2,eind3,eind4];
    if strcmp(eind,'')
        eind = 'N/A';
    end
    evpow = app.EIpower.Value;


    % From Partial Correlation Coefficients
    rpow = app.PCCpower.Value;
    
    % From threshold procession
    thact = [];
    thres1 = [];
    thres2 = [];
    thres3 = [];   

    % From Polynomial Fit
    polyfit = app.PolyfitButton.Value;
    
    % From Saving
    opt = logical(app.OPTButton.Value);
    rnk = logical(app.RNKButton.Value);
    raw = logical(app.RAWButton.Value);
    fc = logical(app.FCButton.Value);
    
    % From Toolbar
    time = logical(app.time.State);
        
    % From Parallel Pool
    ncores = str2double(app.ncoresbox.Value);
    
    % Saving Job information
    job(nj).name = name;
    job(nj).ffts = ffts;
    job(nj).sr = min(cat(1,ffts.sr));
    job(nj).maxf = maxf;
    job(nj).minf = minf;
    job(nj).ftarg = ftarg;
    job(nj).osc = osc;

    job(nj).exclude = excfreq;    
    job(nj).singf = singf;
    job(nj).phsing = phsing;
    job(nj).mult = mult;
    job(nj).phmult = phmult;    

    if time
        job(nj).schedule.state = true;
        job(nj).schedule.days = schedule.days;
        job(nj).schedule.tstart = schedule.tstart;
        job(nj).schedule.tend = schedule.tend;
    else
        job(nj).schedule.state = false;
    end

    job(nj).bivar = bivar;
    job(nj).ev = ev;
    job(nj).model = model;
    job(nj).reg = reg;
    job(nj).outlier = outlier;
    job(nj).remmethod = remmethod;

    job(nj).zrank = zrank;
    job(nj).trank = trank;
    job(nj).ranking = ranking;

    job(nj).eind = eind;
    job(nj).evpow = evpow;
    job(nj).rpow = rpow;
    
    job(nj).thres.active = thact;
    job(nj).thres.thres1 = thres1;
    job(nj).thres.thres2 = thres2;
    job(nj).thres.thres3 = thres3;

    job(nj).polyfit = polyfit;
    
    job(nj).opt = opt;
    job(nj).rnk = rnk;
    job(nj).raw = raw;
    job(nj).fc = fc;
    job(nj).ncores = ncores;
    
    app.job = job;
    
end