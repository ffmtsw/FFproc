function final = FFproc_TFwmedian(job,target_freq,ntseg,tseg,act_chan,time,src,RNK,RAW,EV,TF,FC)
    
    % Number of time series per job
    nts = numel(job.ffts);
    % Maximum of ammount of data to be considered for estimations (1-100%)
    switch src
        case 'SLD'
            ranking = job.ranking;
        case 'ALL'
            ranking = 100;
    end
    
    % Magnetic reference channels
    magref = zeros(numel(job.ffts),2);
    for j = 1:numel(job.ffts)
        magref(j,1) = job.ffts(j).rsbx;
        magref(j,2) = job.ffts(j).rsby;
    end    
        
    % Index structure for fitted estimations
    id = struct;    
    cf = struct;
    % Create FFsave structure to pre-allocate fields.
    FFsave = create_ffsave;

    % Calculating the best transfer functions for each target frequency
    % based on a polynomial fitting, for each time series (nts)
    for a = 1:nts        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                 IMPEDANCES ESTIMATION                     %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Deleting previous estimations and preallocating variables
        clear Zxx Zxy Zyx Zyy Zxx_Err Zxy_Err Zyx_Err Zyy_Err  
        % If Ex and Ey are active, Impedances and Phase Tensor components
        % can be calculated
        if act_chan{a}(1) && act_chan{a}(2)
            Zxx = RAW(ranking).Zxx(:,a);
            Zxy = RAW(ranking).Zxy(:,a);
            Zyx = RAW(ranking).Zyx(:,a);
            Zyy = RAW(ranking).Zyy(:,a);
            Zxx_Err = RAW(ranking).Zxx_Err(:,a);
            Zxy_Err = RAW(ranking).Zxy_Err(:,a);
            Zyx_Err = RAW(ranking).Zyx_Err(:,a);
            Zyy_Err = RAW(ranking).Zyy_Err(:,a);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%              PHASE TENSOR ESTIMATION                  %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            phi11 = RAW(ranking).phi11(:,a);      phi12 = RAW(ranking).phi12(:,a);
            phi21 = RAW(ranking).phi21(:,a);      phi22 = RAW(ranking).phi22(:,a);
            beta = RAW(ranking).beta(:,a);        beta_Err = RAW(ranking).beta_Err(:,a);
            alpha = RAW(ranking).alpha(:,a);      alpha_Err = RAW(ranking).alpha_Err(:,a);
            theta = RAW(ranking).theta(:,a);      theta_Err = RAW(ranking).theta_Err(:,a);
            lambda = RAW(ranking).lambda(:,a);    lambda_Err = RAW(ranking).lambda_Err(:,a);
            phimax = RAW(ranking).phimax(:,a);    phimax_Err = RAW(ranking).phimax_Err(:,a);
            phimin = RAW(ranking).phimin(:,a);    phimin_Err = RAW(ranking).phimin_Err(:,a);                
        else
            % If impedances can not be calculated, fill Z components with 
            % NaN
            Zxx = NaN(numel(target_freq),1);            Zxx_Err = NaN(numel(target_freq),1);
            Zxy = NaN(numel(target_freq),1);            Zxy_Err = NaN(numel(target_freq),1);
            Zyx = NaN(numel(target_freq),1);            Zyx_Err = NaN(numel(target_freq),1);
            Zyy = NaN(numel(target_freq),1);            Zyy_Err = NaN(numel(target_freq),1);
            % If impedances can not be calculated, fill PT components with 
            % NaN
            phi11 = NaN(numel(target_freq),1);          phi12 = NaN(numel(target_freq),1);
            phi21 = NaN(numel(target_freq),1);          phi22 = NaN(numel(target_freq),1);
            beta = NaN(numel(target_freq),1);           beta_Err = NaN(numel(target_freq),1);
            alpha = NaN(numel(target_freq),1);          alpha_Err = NaN(numel(target_freq),1);
            theta = NaN(numel(target_freq),1);          theta_Err = NaN(numel(target_freq),1);
            lambda = NaN(numel(target_freq),1);         lambda_Err = NaN(numel(target_freq),1);
            phimax = NaN(numel(target_freq),1);         phimax_Err = NaN(numel(target_freq),1);
            phimin = NaN(numel(target_freq),1);         phimin_Err = NaN(numel(target_freq),1); 
        end     
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                   TIPPER ESTIMATION                       %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Deleting previous estimations and preallocating variables
        clear txz tyz txz_Err tyz_Err        
        % If Bz is active, Tipper components can be calculated
        if act_chan{a}(5)
            % Extracting info from Results structure variable
            txz = RAW(ranking).txz(:,a);
            tyz = RAW(ranking).tyz(:,a);
            txz_Err = RAW(ranking).txz_Err(:,a);
            tyz_Err = RAW(ranking).tyz_Err(:,a);
        else
            % If tippers can not be calculated, fill T components with NaN
            txz = NaN(numel(target_freq),1);            txz_Err = NaN(numel(target_freq),1);
            tyz = NaN(numel(target_freq),1);            tyz_Err = NaN(numel(target_freq),1);
        end       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                 SAVING FINAL ESTIMATIONS                  %%%%
        %%%%                   FOR EACH TIME SERIES                    %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        % Saving number of frequencies, frequency and period vector
        FFsave(a).nfreq = numel(target_freq);
        FFsave(a).freq = target_freq;
        FFsave(a).per = 1./(target_freq);
        % Saving transfer functions in FFsave structure: Z
        FFsave(a).Zxx = Zxx;                            FFsave(a).Zxx_Err = Zxx_Err;
        FFsave(a).Zxy = Zxy;                            FFsave(a).Zxy_Err = Zxy_Err;
        FFsave(a).Zyx = Zyx;                            FFsave(a).Zyx_Err = Zyx_Err;
        FFsave(a).Zyy = Zyy;                            FFsave(a).Zyy_Err = Zyy_Err;
        % Saving transfer functions in FFsave structure: T
        FFsave(a).txz = txz;                            FFsave(a).txz_Err = txz_Err;
        FFsave(a).tyz = tyz;                            FFsave(a).tyz_Err = tyz_Err;
        % Saving transfer functions in results structure: PT
        FFsave(a).phi11 = phi11;                        FFsave(a).phi12 = phi12;
        FFsave(a).phi21 = phi21;                        FFsave(a).phi22 = phi22;
        FFsave(a).beta = beta;                          FFsave(a).beta_Err = beta_Err;
        FFsave(a).alpha = alpha;                        FFsave(a).alpha_Err = alpha_Err;
        FFsave(a).theta = theta;                        FFsave(a).theta_Err = theta_Err;
        FFsave(a).lambda = lambda;                      FFsave(a).lambda_Err = lambda_Err;
        FFsave(a).phimax = phimax;                      FFsave(a).phimax_Err = phimax_Err;
        FFsave(a).phimin = phimin;                      FFsave(a).phimin_Err = phimin_Err;       
    end
        
    % Writing extra info to the resultant FFsave variable:
    % Name, coordinates (lonlat, UTM) and elevation.
    FFsave = FFproc_OPT_metainfo(FFsave,job);
    
    % Resistivity and phase calculations, escaling and corrections
    FFsave = FFproc_OPT_rhoandphi(FFsave,job);   
    
    % Saving Eigenvalues information
    try
        EV = FFproc_OPT_saveEV(job,target_freq,ntseg,tseg,id,EV); 
    catch
        EV = [];
    end
    
    % Saving in final structure
    final.job = job;                    % job parameters
    final.target_freq = target_freq;    % Evaluated target frequencies
    final.act_chan = act_chan;          % Active channels for processing
    final.time = time;                  % Date and time for each time segment
    final.FFsave = FFsave;              % Structure with MT-responses
    final.RNK = RNK;                    % Ranking parameters with raw TF
    final.RAW = RAW;                    % Estimations for all percentages
    final.EV = EV;                      % Structure with Eigenvalues
    final.id = id;                      % %selected windows for each tf
    final.cf = cf;                      % Polyomial fitting
    final.TF = TF;                      % Transfer functions
    final.FC = FC;                      % Fourier Coeficients    

end