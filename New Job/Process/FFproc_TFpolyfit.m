function final = FFproc_TFpolyfit(job,target_freq,ntseg,tseg,act_chan,time,RNK,RAW,EV,TF,FC)
    
    % Number of time series per job
    nts = numel(job.ffts);
    % Maximum of ammount of data to be considered for estimations (1-100%)
    ranking = job.ranking;
    
    % Magnetic reference channels
    magref = zeros(numel(job.ffts),2);
    for j = 1:numel(job.ffts)
        magref(j,1) = job.ffts(j).rsbx;
        magref(j,2) = job.ffts(j).rsby;
    end    
    
    try
        % Polynomial order defined by FFproc_RAW
        polyord = job.polyord;
    catch
        % Polynomial order (the polynom order will be 1.25 times of number
        % of decades along target frequencies range)
        polyord = round(1.25*(ceil(max(log10(target_freq))) - floor(min(log10(target_freq)))));   
        if polyord > 9
            polyord = 9;
        end
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
        %%%%                   IMPEDANCES ESTIMATION                   %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Deleting previous estimations and preallocating variables
        clear ZXX ZXY ZYX ZYY ZXX_Err ZXY_Err ZYX_Err ZYY_Err
        clear Zxx Zxy Zyx Zyy Zxx_Err Zxy_Err Zyx_Err Zyy_Err  
        ZXX = zeros(numel(target_freq),ranking);        ZXX_Err = zeros(numel(target_freq),ranking);
        ZXY = zeros(numel(target_freq),ranking);        ZXY_Err = zeros(numel(target_freq),ranking);
        ZYX = zeros(numel(target_freq),ranking);        ZYX_Err = zeros(numel(target_freq),ranking);
        ZYY = zeros(numel(target_freq),ranking);        ZYY_Err = zeros(numel(target_freq),ranking);
        % If Ex and Ey are active, Impedances and Phase Tensor components
        % can be calculated
        if act_chan{a}(1) && act_chan{a}(2)
            for f = 1:numel(target_freq)
                for r = 1:ranking
                    ZXX(f,r) = RAW(r).Zxx(f,a);
                    ZXY(f,r) = RAW(r).Zxy(f,a);
                    ZYX(f,r) = RAW(r).Zyx(f,a);
                    ZYY(f,r) = RAW(r).Zyy(f,a);
                    ZXX_Err(f,r) = RAW(r).Zxx_Err(f,a);
                    ZXY_Err(f,r) = RAW(r).Zxy_Err(f,a);
                    ZYX_Err(f,r) = RAW(r).Zyx_Err(f,a);
                    ZYY_Err(f,r) = RAW(r).Zyy_Err(f,a);                    
                end
            end
            % Polynom fitting for Z
            [Zxx,Zxx_Err,id(a).Zxx_Re,id(a).Zxx_Im,cf(a).Zxx_Re,cf(a).Zxx_Im] = polynomialfit(target_freq,ZXX,ZXX_Err,'ZXX',polyord);
            [Zxy,Zxy_Err,id(a).Zxy_Re,id(a).Zxy_Im,cf(a).Zxy_Re,cf(a).Zxy_Im] = polynomialfit(target_freq,ZXY,ZXY_Err,'ZXY',polyord);
            [Zyx,Zyx_Err,id(a).Zyx_Re,id(a).Zyx_Im,cf(a).Zyx_Re,cf(a).Zyx_Im] = polynomialfit(target_freq,ZYX,ZYX_Err,'ZYX',polyord);
            [Zyy,Zyy_Err,id(a).Zyy_Re,id(a).Zyy_Im,cf(a).Zyy_Re,cf(a).Zyy_Im] = polynomialfit(target_freq,ZYY,ZYY_Err,'ZYY',polyord);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%                PHASE TENSOR ESTIMATION                %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Deleting previous estimations and preallocating variables
            clear Z Z_Err PT PT_Err
            phi11 = zeros(numel(target_freq),1);        phi12 = zeros(numel(target_freq),1);
            phi21 = zeros(numel(target_freq),1);        phi22 = zeros(numel(target_freq),1);
            alpha = zeros(numel(target_freq),1);        beta = zeros(numel(target_freq),1);
            theta = zeros(numel(target_freq),1);        lambda = zeros(numel(target_freq),1);
            phimax = zeros(numel(target_freq),1);       phimin = zeros(numel(target_freq),1);
            sig_min = zeros(numel(target_freq),1);      sig_max = zeros(numel(target_freq),1);
            sig_alpha = zeros(numel(target_freq),1);    sig_beta = zeros(numel(target_freq),1);
            sig_theta = zeros(numel(target_freq),1);    sig_lambda = zeros(numel(target_freq),1);            
            % Phase Tensor calculations for each target frequency
            for f = 1:numel(target_freq)
                Z = [Zxx(f);Zyx(f);Zxy(f);Zyy(f)];
                Z_err = [Zxx_Err(f);Zyx_Err(f);Zxy_Err(f);Zyy_Err(f)];
                PT = FFproc_PT(Z);
                PT_err = FFproc_PT_err(Z,Z_err);
                phi11(f,1) = PT.phi11;                  phi12(f,1) = PT.phi12;
                phi21(f,1) = PT.phi21;                  phi22(f,1) = PT.phi22;
                alpha(f,1) = PT.alpha;                  beta(f,1) = PT.beta;
                theta(f,1) = PT.theta;                  lambda(f,1) = PT.lambda;
                phimax(f,1) = PT.phimax;                phimin(f,1) = PT.phimin;
                sig_min(f,1) = PT_err(1);               sig_max(f,1) = PT_err(2);
                sig_alpha(f,1) = PT_err(4);             sig_beta(f,1) = PT_err(3);
                sig_theta(f,1) = PT_err(5);             sig_lambda(f,1) = PT_err(6);
            end
        else
            % If impedances can not be calculated, fill Z components with 
            % NaN
            Zxx = NaN(numel(target_freq),1);            Zxx_Err = NaN(numel(target_freq),1);
            Zxy = NaN(numel(target_freq),1);            Zxy_Err = NaN(numel(target_freq),1);
            Zyx = NaN(numel(target_freq),1);            Zyx_Err = NaN(numel(target_freq),1);
            Zyy = NaN(numel(target_freq),1);            Zyy_Err = NaN(numel(target_freq),1);
            id(a).Zxy_Re = 2*ones(numel(target_freq),1);
            id(a).Zxy_Im = 2*ones(numel(target_freq),1);  
            id(a).Zyx_Re = 2*ones(numel(target_freq),1);
            id(a).Zyy_Im = 2*ones(numel(target_freq),1);
            cf(a).Zxy_Re = [];
            cf(a).Zxy_Im = [];
            cf(a).Zyx_Re = [];
            cf(a).Zyx_Im = [];
            % If impedances can not be calculated, fill PT components with 
            % NaN
            phi11 = NaN(numel(target_freq),1);          phi12 = NaN(numel(target_freq),1);
            phi21 = NaN(numel(target_freq),1);          phi22 = NaN(numel(target_freq),1);
            alpha = NaN(numel(target_freq),1);          beta = NaN(numel(target_freq),1);
            theta = NaN(numel(target_freq),1);          lambda = NaN(numel(target_freq),1);
            phimax = NaN(numel(target_freq),1);         phimin = NaN(numel(target_freq),1);
            sig_min = NaN(numel(target_freq),1);        sig_max = NaN(numel(target_freq),1);
            sig_alpha = NaN(numel(target_freq),1);      sig_beta = NaN(numel(target_freq),1);
            sig_theta = NaN(numel(target_freq),1);      sig_lambda = NaN(numel(target_freq),1);
        end     
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                     TIPPER ESTIMATION                     %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Deleting previous estimations and preallocating variables
        clear TXZ TYZ TXZ_Err TYZ_Err  
        clear txz tyz txz_Err tyz_Err        
        TXZ = zeros(numel(target_freq),ranking);        TXZ_Err = zeros(numel(target_freq),ranking);
        TYZ = zeros(numel(target_freq),ranking);        TYZ_Err = zeros(numel(target_freq),ranking);  
        % If Bz is active, Tipper components can be calculated
        if act_chan{a}(5)
            % Extracting info from Results structure variable
            for f = 1:numel(target_freq)  % Number of time segments per job
                for r = 1:ranking  % Eigenvalue percentage
                    TXZ(f,r) = RAW(r).txz(f,a);
                    TYZ(f,r) = RAW(r).tyz(f,a);
                    TXZ_Err(f,r) = RAW(r).txz_Err(f,a);
                    TYZ_Err(f,r) = RAW(r).tyz_Err(f,a);
                end
            end
            % Polynom fitting for T
            [txz,txz_Err,id(a).txz_Re,id(a).txz_Im,cf(a).txz_Re,cf(a).txz_Im] = polynomialfit(target_freq,TXZ,TXZ_Err,'TXZ',polyord);
            [tyz,tyz_Err,id(a).tyz_Re,id(a).tyz_Im,cf(a).tyz_Re,cf(a).tyz_Im] = polynomialfit(target_freq,TYZ,TYZ_Err,'TYZ',polyord);
        else
            % If tippers can not be calculated, fill T components with NaN
            txz = NaN(numel(target_freq),1);            txz_Err = NaN(numel(target_freq),1);
            tyz = NaN(numel(target_freq),1);            tyz_Err = NaN(numel(target_freq),1);
            id(a).txz_Re = 2*ones(numel(target_freq),1);   
            id(a).txz_Im = 2*ones(numel(target_freq),1);   
            id(a).tyz_Re = 2*ones(numel(target_freq),1);
            id(a).tyz_Im = 2*ones(numel(target_freq),1);            
            cf(a).txz_Re = [];
            cf(a).txz_Im = [];
            cf(a).tyz_Re = [];
            cf(a).tyz_Im = [];
        end       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%                   SAVING FINAL ESTIMATIONS                %%%%
        %%%%                     FOR EACH TIME SERIES                  %%%%
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
        FFsave(a).beta = beta;                          FFsave(a).beta_Err = sig_beta;
        FFsave(a).alpha = alpha;                        FFsave(a).alpha_Err = sig_alpha;
        FFsave(a).theta = theta;                        FFsave(a).theta_Err = sig_theta;
        FFsave(a).lambda = lambda;                      FFsave(a).lambda_Err = sig_lambda;
        FFsave(a).phimax = phimax;                      FFsave(a).phimax_Err = sig_max;
        FFsave(a).phimin = phimin;                      FFsave(a).phimin_Err = sig_min;       
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