function FFsave = ffsave_create

    % Initializing variable
    FFsave = struct;
    % Header
    FFsave.ID = [];
    FFsave.site = [];         % Site name
    FFsave.lonlat = [];       % Geographic position
    FFsave.UTM = [];          % UTM Coordinates
    FFsave.z = [];            % Elevation
    FFsave.nfreq = [];        % Number of frequencies
    FFsave.freq = [];         % Frequencies
    FFsave.per = [];          % Periods
    % Impedances
    FFsave.Zxx = [];        FFsave.Zxx_Err = [];
    FFsave.Zxy = [];        FFsave.Zxy_Err = [];
    FFsave.Zyx = [];        FFsave.Zyx_Err = [];
    FFsave.Zyy = [];        FFsave.Zyy_Err = [];
    % Resistivities
    FFsave.rhoxx = [];      FFsave.rhoxy = [];    
    FFsave.rhoyx = [];      FFsave.rhoyy = [];
    FFsave.rhoxx_Err = [];  FFsave.rhoxy_Err = [];    
    FFsave.rhoyx_Err = [];  FFsave.rhoyy_Err = [];
    % Phases
    FFsave.phixx = [];      FFsave.phixy = [];    
    FFsave.phiyx = [];      FFsave.phiyy = [];
    FFsave.phixx_Err = [];  FFsave.phixy_Err = [];    
    FFsave.phiyx_Err = [];  FFsave.phiyy_Err = [];
    % Tipper
    FFsave.txz = [];        FFsave.txz_Err = [];
    FFsave.tyz = [];        FFsave.tyz_Err = [];
    % Phase Tensor
    FFsave.phi11 = [];    
    FFsave.phi12 = [];
    FFsave.phi21 = [];
    FFsave.phi22 = [];
    FFsave.beta = [];       FFsave.beta_Err = [];
    FFsave.alpha = [];      FFsave.alpha_Err = [];
    FFsave.theta = [];      FFsave.theta_Err = [];
    FFsave.lambda = [];     FFsave.lambda_Err = [];
    FFsave.phimax = [];     FFsave.phimax_Err = [];
    FFsave.phimin = [];     FFsave.phimin_Err = [];
    % Horizontal Magnetic Transfer Functions
    FFsave.txx = [];        FFsave.txy = [];
    FFsave.tyx = [];        FFsave.tyy = [];
    FFsave.t_id = [];
    
end