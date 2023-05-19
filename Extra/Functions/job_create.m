function job = job_create

    % Name of sites included in the job
    job.name = [];
    % Time series structure
    job.ffts = [];
    % Sampling rate (Hz)
    job.sr = [];
    % Maximum frequency range to be estimated
    job.maxf = [];
    % Minimum frequency range to be estimated
    job.minf = [];
    % Target frequencies per decade
    job.ftarg = [];
    % Number of oscillations to calulate time segments
    job.osc = [];

    % Remove spectral lines
    job.exclude = [];  
    % Single frequencies
    job.singf = [];
    % Width of filter for single frequencies
    job.phsing = [];
    % Multiples
    job.mult = [];
    % Width of filter for multiples
    job.phmult = [];

    % Schedule
    job.schedule.state = false;
    % Days
    job.schedule.days = [];
    % Start time
    job.schedule.tstart = [];
    % End time
    job.schedule.tend = [];

    % Transfer function estimation approach
    job.bivar = [];
    job.ev = [];
    % Noise Mode
    job.model = [];
    % Linear regression
    job.reg = [];
    % Outlier dettection
    job.outlier = [];
    % Detection method
    job.remmethod = [];
    % Significance Level
    job.alpha = [];
    
    % Ranking Criterion for Z
    job.zrank = [];
    % Ranking Criterion for T
    job.trank = [];
    % Ranking ammount (%)
    job.ranking = [];

    % Eigenvalue Index
    job.eind = [];
    % Power for Eigenvalue Index
    job.evpow = [];
    % Power of Partial Correlation Coefficients
    job.rpow = [];

    % Threshold Processing
    job.thres.active = [];
    job.thres.thres1 = [];
    job.thres.thres2 = [];
    job.thres.thres3 = [];

    % Polynomial Fit
    job.polyfit = [];

    % OPT file
    job.opt = [];
    % RAW file
    job.raw = [];
    % Fourier Coefficients file
    job.fc = [];
    % Ranking Parameters file
    job.rnk = false;

    % Number of MATLAB Workers
    job.ncores = [];

end