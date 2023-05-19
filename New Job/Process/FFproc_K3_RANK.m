% function FFproc_K3_RNK is the 3/3 core functions of FFproc_Process.m
% This function evaluates the best time segments based on the Eigenvalue
% criterion, calculating an Eigenvalue Index (Hering, 2015; Hering 2019)
% and Partial Coherences.
% It calculates transfer functions (Z, T and PT) for those best time 
% segments. 
% 
% A weighted median calculation is used to obtain a robust
% estimation of each Z and T components. Errors are calculated based on 
% Häuserer(2010).
% version 1.0 / ?????2019 / ph
% version 1.1 / ?????2019 / at
% version 2.0 / 12apr2020 / cc   Modified to be implemented within the new
%                                version of FFMT (2.0). Name has changed
%                                from TF_aut.
% version 2.1 / 23nov2020 / cc   Modified to take into account active
%                                channels (for time series without electric
%                                and magnetic Bz channels, e.g. RR stations
%                                only recording Bx and By channels. If
%                                channels Ex, Ey and/or Bz are not active,
%                                the result is replaced by zeros
% version 3.0 / 09aug2022 / cc   Modified to take into account the partial
%                                correlation coefficients calculated in K1.
%                                If the PCC option is active, the ranking
%                                of the time segments is based on the
%                                magnitude of the PCC for each Impedance
%                                and Tipper components.

function [results,EigVal,RNK] = FFproc_K3_RANK(job,ntseg,act_chan,TF,RNK,parproc)

    warning('off','all')
    
    % Number of time series per job
    nts = numel(job.ffts);

    % Bivariate estimation of Z and T when using single site
    bivar = job.bivar;

    % Number of workers
    if parproc         
        ncores = job.ncores;
    else
        ncores = 0;
    end

    % Eigenvalue Criterion
    EVC = job.eind;
    
    % Power of EI
    evpow = job.evpow;

    % Power of PCC
    rpow = job.rpow;

    % Outlier detection and removal
    tfout = job.outlier;    

    % If Mahalanobis outlier method is chosen, no outlier detection is
    % performed during this kernel. Only for testing (24mar2023 / cc)
    remmethod = job.remmethod;
    switch remmethod
        case 'maha'
            tfout = false;
    end

    % Magnetic reference channels
    magref = zeros(numel(job.ffts),2);
    for i = 1:numel(job.ffts)
        magref(i,1) = job.ffts(i).rsbx;
        magref(i,2) = job.ffts(i).rsby;
    end

    % ------------------------------------------------------------------- %
    %                   Partial coherences extraction                     %
    % ------------------------------------------------------------------- %        
    % Extracting partial coherences computed in K2
    Zxx_pcorr = cat(1,TF.rZxx);
    Zxy_pcorr = cat(1,TF.rZxy);
    Zyx_pcorr = cat(1,TF.rZyx);    
    Zyy_pcorr = cat(1,TF.rZyy);
    txz_pcorr = cat(1,TF.rtxz);
    tyz_pcorr = cat(1,TF.rtyz);
           
    % ------------------------------------------------------------------- %
    %                       Eigenvalues extraction                        %
    % ------------------------------------------------------------------- %        
    % Extracting eigenvalues computed in K2    
    if bivar
        % If bivariate approach is choosen, only partial coherences are
        % taken into account as weights. Therefore, Eigenvalue Indices are
        % ommited and defined as 1.
        EI0 = zeros(ntseg,1);
        EI1 = zeros(ntseg,1);
        EI2 = zeros(ntseg,1);
        EI3 = zeros(ntseg,1);
        EI4 = zeros(ntseg,1);
        EI = ones(ntseg,1);
    else
        % Extracting Eigenvalues for all time segments (ntseg x n_act_chan)
        EVmat = cat(2,TF.EVal).';
       
        % Eigenvalue 0: |λ2|/|λ1|
        if any(ismember(string(EVC),'EI0')) || any(ismember(string(EVC),'EV0'))
            EI0 = calcevindex(EVmat,0,tfout,remmethod);
        else
            EI0 = ones(ntseg,1);
        end        
        % Eigenvalue 1: 1 - |λ3|/|λ2|
        if any(ismember(string(EVC),'EI1')) || any(ismember(string(EVC),'EV1'))
            EI1 = calcevindex(EVmat,1,tfout,remmethod);
        else
            EI1 = ones(ntseg,1);
        end
        % Eigenvalue 2: 1/|λ3|
        if any(ismember(string(EVC),'EI2')) || any(ismember(string(EVC),'EV2'))
            EI2 = calcevindex(EVmat,2,tfout,remmethod);
        else
            EI2 = ones(ntseg,1);
        end
        % Eigenvalue 3: |λ2|/|λ3| - |λ3|^2
        if any(ismember(string(EVC),'EI3')) || any(ismember(string(EVC),'EV3'))
            EI3 = calcevindex(EVmat,3,tfout,remmethod);
        else
            EI3 = ones(ntseg,1);
        end
        % Eigenvalue 4: |λ2|/|λ3| - |λ1|/|λ2|
        if any(ismember(string(EVC),'EI4')) || any(ismember(string(EVC),'EV4'))
            EI4 = calcevindex(EVmat,4,tfout,remmethod);
        else
            EI4 = ones(ntseg,1);
        end
        % Calculating Eigenvalue Index (EI) as Hering et al. (2023) for each 
        % time window (w)
        EI = (EI0.^evpow).*(EI1.^evpow).*(EI2.^evpow).*(EI3.^evpow).*(EI4.^evpow);
    end  

    % ------------------------------------------------------------------- %
    %                                                                     %
    %                 Saving criterion/ranking parameters                 %
    %                                                                     %
    %         This part will be only executed if the user chooses         %
    %             to save the RNK file. This will contain RAW             %
    %     parameters that can be used aferwards during post-processing    %
    %                   using FFproc - RNK application.                   %
    %          The user, based on visual inspection, will be able         %
    %           to change specific parameters (EI, PCOH, others)          %
    %                which yields in better estimations of TF             %
    %                                                                     %
    % ------------------------------------------------------------------- %  
    if job.rnk
        for w = 1:ntseg
            % Transfer functions
            RNK(w).Zxx = TF(w).Zxx;
            RNK(w).Zxy = TF(w).Zxy;
            RNK(w).Zyx = TF(w).Zyx;
            RNK(w).Zyy = TF(w).Zyy;
            RNK(w).txz = TF(w).txz;
            RNK(w).tyz = TF(w).tyz;
            % Partial Coherences for Impedances and Tipper
            RNK(w).Zxx_pcorr = Zxx_pcorr(w,:);
            RNK(w).Zxy_pcorr = Zxy_pcorr(w,:);
            RNK(w).Zyx_pcorr = Zyx_pcorr(w,:);
            RNK(w).Zyy_pcorr = Zyy_pcorr(w,:);
            RNK(w).txz_pcorr = txz_pcorr(w,:);
            RNK(w).tyz_pcorr = tyz_pcorr(w,:);  
            % Saving Eigenvalues
            RNK(w).EVal = TF(w).EVal;
            % Eigenvalue Indices
            RNK(w).EI0 = EI0(w);
            RNK(w).EI1 = EI1(w);
            RNK(w).EI2 = EI2(w);
            RNK(w).EI3 = EI3(w);
            RNK(w).EI4 = EI4(w);
            RNK(w).EI = EI(w);           
        end
    else
        RNK = struct;
    end

    % ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| %
    % VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV %
    %                                                                     %
    %                    Start of time segment ranking                    %
    %                                                                     %
    % ------------------------------------------------------------------- % 

    % Pre-allocating space for Impedance Z
    EIPC_Zxx = NaN(size(Zxx_pcorr));    EIPC_Zxx_ind = NaN(size(Zxx_pcorr));
    EIPC_Zxy = NaN(size(Zxy_pcorr));    EIPC_Zxy_ind = NaN(size(Zxy_pcorr));
    EIPC_Zyx = NaN(size(Zyx_pcorr));    EIPC_Zyx_ind = NaN(size(Zyx_pcorr));
    EIPC_Zyy = NaN(size(Zyy_pcorr));    EIPC_Zyy_ind = NaN(size(Zyy_pcorr));
    % Pre-allocating space for Tipper T
    EIPC_txz = NaN(size(txz_pcorr));    EIPC_txz_ind = NaN(size(txz_pcorr));
    EIPC_tyz = NaN(size(tyz_pcorr));    EIPC_tyz_ind = NaN(size(tyz_pcorr));

    % Option to choose which approach to be used:
    % To take into account the EigenValue Index, EV should be true,
    % otherwise the values of EI will be taken as 1.
    % To take into account the PCC for each component of Z and T, PCC 
    % should be true, otherwise the PCC will be taken as 1
    switch job.zrank
        case 'EV'
            EI_Z = EI;
            rpow_Z = 0;
        case 'PCC'
            EI_Z = EI.^0;
            rpow_Z = job.rpow;
        case 'EVPCC'
            EI_Z = EI;
            rpow_Z = job.rpow;
    end
    switch job.trank
        case 'EV'
            EI_T = EI;
            rpow_T = 0;
        case 'PCC'
            EI_T = EI.^0;
            rpow_T = job.rpow;
        case 'EVPCC'
            EI_T = EI;
            rpow_T = job.rpow;
    end

    % Scaling Eigenvalue index by the partial correlation coefficients used
    % to calculate each component of the Impedance Tensor and Tipper
    % Vectors. EI is the same for all the components during the same time
    % segment, but Partial coherences are not. Therefore, different time
    % windows may be taken to calculate Impedance and Tipper components
    for a = 1:nts  
        % ZXX
        EP_Zxx = EI_Z.*Zxx_pcorr(:,a).^rpow_Z;
        [EIPC_Zxx(:,a),EIPC_Zxx_ind(:,a)] = sort(EP_Zxx,'descend');        
        % ZXY
        EP_Zxy = EI_Z.*Zxy_pcorr(:,a).^rpow_Z;
        [EIPC_Zxy(:,a),EIPC_Zxy_ind(:,a)] = sort(EP_Zxy,'descend');          
        % ZYX
        EP_Zyx = EI_Z.*Zyx_pcorr(:,a).^rpow_Z;
        [EIPC_Zyx(:,a),EIPC_Zyx_ind(:,a)] = sort(EP_Zyx,'descend');        
        % ZYY
        EP_Zyy = EI_Z.*Zyy_pcorr(:,a).^rpow_Z;
        [EIPC_Zyy(:,a),EIPC_Zyy_ind(:,a)] = sort(EP_Zyy,'descend');    
        % TXZ
        EP_txz = EI_T.*txz_pcorr(:,a).^rpow_T;
        [EIPC_txz(:,a),EIPC_txz_ind(:,a)] = sort(EP_txz,'descend');    
        % TYZ
        EP_tyz = EI_T.*tyz_pcorr(:,a).^rpow_T;
        [EIPC_tyz(:,a),EIPC_tyz_ind(:,a)] = sort(EP_tyz,'descend');  
    end

    % ------------------------------------------------------------------- %
    %                  Sorting Transfer Functions and EV                  %
    % ------------------------------------------------------------------- % 
    % Preallocating space for sorted TF and Eigenvalues
    ZXX = cat(1,TF.Zxx);    Zxx = complex(NaN(ntseg,nts),NaN(ntseg,nts));     Zxx_EV = cell(1,nts);
    ZXY = cat(1,TF.Zxy);    Zxy = complex(NaN(ntseg,nts),NaN(ntseg,nts));     Zxy_EV = cell(1,nts);
    ZYX = cat(1,TF.Zyx);    Zyx = complex(NaN(ntseg,nts),NaN(ntseg,nts));     Zyx_EV = cell(1,nts);
    ZYY = cat(1,TF.Zyy);    Zyy = complex(NaN(ntseg,nts),NaN(ntseg,nts));     Zyy_EV = cell(1,nts);
    TXZ = cat(1,TF.txz);    Txz = complex(NaN(ntseg,nts),NaN(ntseg,nts));     txz_EV = cell(1,nts);
    TYZ = cat(1,TF.tyz);    Tyz = complex(NaN(ntseg,nts),NaN(ntseg,nts));     tyz_EV = cell(1,nts);
    
    % From now on, the Transfer Functions Z and T  and EV are already 
    % sorted based on the Eigenvalue-Partial Coherence (EIPC) Index. 
    % Each station, each component of TF and EV are sorted.
    for a = 1:nts
        % Sorting transfer functions: Impedance
        Zxx(:,a) = ZXX(EIPC_Zxx_ind(:,a),a);
        Zxy(:,a) = ZXY(EIPC_Zxy_ind(:,a),a);
        Zyx(:,a) = ZYX(EIPC_Zyx_ind(:,a),a);
        Zyy(:,a) = ZYY(EIPC_Zyy_ind(:,a),a);
        % Sorting transfer functions: Tipper
        Txz(:,a) = TXZ(EIPC_txz_ind(:,a),a);
        Tyz(:,a) = TYZ(EIPC_tyz_ind(:,a),a);
        if ~bivar
            % Sorting EigenValues: Impedance
            Zxx_EV{a} = cat(2,TF(EIPC_Zxx_ind(:,a)).EVal).';
            Zxy_EV{a} = cat(2,TF(EIPC_Zxy_ind(:,a)).EVal).';
            Zyx_EV{a} = cat(2,TF(EIPC_Zyx_ind(:,a)).EVal).';
            Zyy_EV{a} = cat(2,TF(EIPC_Zyy_ind(:,a)).EVal).';
            % Sorting EigenValues: Tipper
            txz_EV{a} = cat(2,TF(EIPC_txz_ind(:,a)).EVal).';
            tyz_EV{a} = cat(2,TF(EIPC_tyz_ind(:,a)).EVal).';
        end
    end
    % Release memory of non sorted transfer functions
    clear ZXX ZXY ZYX ZYY TXZ TYZ

    % Allocating structure 
    results = struct;
    segments = struct;

    % Estimation for different ammount of data (1 - 100%)
    parfor (r = 1:100,ncores)
%     for r = 1:100

        % ntseg:    number of total time segments
        % seg_eval: number of time segments evaluated for each ranking
        seg_eval = ceil(r*ntseg/100);
        if seg_eval > ntseg
            seg_eval = ntseg;            
        end

        % Saving segments used to estimate each percentage
        segments(r).seg_eval = seg_eval;

        %                     ** IMPORTANT NOTE ** 
        % The developers acknowledge that this process is not the fastest,
        % as loops over stations are carried for weights, eigenvalues and 
        % transfer functions separately.
        % However it is easy for the user to follow-up the code and to 
        % understand how the ranking is performed. We calculate rankings as 
        % the following:

        % TF/EV - 1%
        % TF/EV -- 2%
        % TF/EV --- 3%
        % TF/EV ---- 4%
        % TF/EV ----- 5%
        % TF/EV ------------------ r%
        % TF/EV ------------------------------------ user defined ranking %

        % Therefore, we will have as many estimtions for TF and EV as the
        % user defined ranking percentage.
        
        % --------------------------------------------------------------- %
        %                             WEIGHTS                             %
        % --------------------------------------------------------------- % 
        % Preallocating space for weights for each component and all
        % stations
        Wxx = cell(1,nts);
        Wxy = cell(1,nts);
        Wyx = cell(1,nts);
        Wyy = cell(1,nts);
        Wxz = cell(1,nts);
        Wyz = cell(1,nts);  

        for a = 1:nts
            % Extracting weights for each station for each level of ranking 
            % percentage. During this stage, the weights are transformed
            % into a right-sided cosine using a half-tukey window, and
            % giving a weight of 1 for the first 5% of data
            % Developers are trying to define which criteron and if real 
            % weights are better for time window selection. Everything is
            % prepared to include independent weights as soon as this issue
            % is solved.

            Wxx{a} = EIPC_Zxx(1:seg_eval,a); %#ok<*PFBNS> 
            Wxy{a} = EIPC_Zxy(1:seg_eval,a);
            Wyx{a} = EIPC_Zyx(1:seg_eval,a);
            Wyy{a} = EIPC_Zyy(1:seg_eval,a);
            Wxz{a} = EIPC_txz(1:seg_eval,a);
            Wyz{a} = EIPC_tyz(1:seg_eval,a);

            % Instead, create weights by a falf-tukey window, with the
            % first 5% of the window weighted as 1. The tukey window will
            % have the same dimensions as the time segment.
%             Wxx{a} = createweights(seg_eval,seg_eval,'chi',5,1,1);
%             Wxy{a} = Wxx{a};
%             Wyx{a} = Wxx{a};
%             Wyy{a} = Wxx{a};
%             Wxz{a} = Wxx{a};
%             Wyz{a} = Wxx{a};
        end

        % --------------------------------------------------------------- %
        %                           EIGENVALUES                           %
        % --------------------------------------------------------------- % 
        for a = 1:nts
            % If bivariate approach is choosen, Eigenvalues cannot be
            % calculated and therefore we use NaN
            if bivar                              
                Zxx_EVal = complex(NaN,NaN);                 
                Zxy_EVal = complex(NaN,NaN);                  
                Zyx_EVal = complex(NaN,NaN);                
                Zyy_EVal = complex(NaN,NaN);                             
                txz_EVal = complex(NaN,NaN);                   
                tyz_EVal = complex(NaN,NaN);
            else    
                % Ammount of the sime segments chosen to be used for
                % calculations (Eigenvalues)
                Zxx_EVeval = Zxx_EV{a}(1 : seg_eval,:);
                Zxy_EVeval = Zxy_EV{a}(1 : seg_eval,:);
                Zyx_EVeval = Zyx_EV{a}(1 : seg_eval,:);
                Zyy_EVeval = Zyy_EV{a}(1 : seg_eval,:);    
                txz_EVeval = txz_EV{a}(1 : seg_eval,:);
                tyz_EVeval = tyz_EV{a}(1 : seg_eval,:);
    
                % Calculating Weighted Medians for EigenValues for Z 
                if act_chan{a}(1) && act_chan{a}(2)  
                    Zxx_EVal = wmedian(Zxx_EVeval,Wxx{a}).';                 
                    Zxy_EVal = wmedian(Zxy_EVeval,Wxy{a}).';                  
                    Zyx_EVal = wmedian(Zyx_EVeval,Wyx{a}).';                
                    Zyy_EVal = wmedian(Zyy_EVeval,Wyy{a}).';
                else
                    Zxx_EVal = complex(NaN,NaN);                 
                    Zxy_EVal = complex(NaN,NaN);                  
                    Zyx_EVal = complex(NaN,NaN);                
                    Zyy_EVal = complex(NaN,NaN);
                end
                % Calculating Weighted Medians for EigenValues for T
                if act_chan{a}(5)                                          
                    txz_EVal = wmedian(txz_EVeval,Wxz{a}).';                   
                    tyz_EVal = wmedian(tyz_EVeval,Wyz{a}).';
                else
                    txz_EVal = complex(NaN,NaN);                   
                    tyz_EVal = complex(NaN,NaN);
                end
            end
            % Saving Weighted EV for each component
            segments(r).Zxx_EVal(a).data = Zxx_EVal;
            segments(r).Zxy_EVal(a).data = Zxy_EVal;
            segments(r).Zyx_EVal(a).data = Zyx_EVal;
            segments(r).Zyy_EVal(a).data = Zyy_EVal;
            segments(r).txz_EVal(a).data = txz_EVal;
            segments(r).tyz_EVal(a).data = tyz_EVal;
        end
        
        % --------------------------------------------------------------- %
        %                        TRANSFER FUNCTIONS                       %
        % --------------------------------------------------------------- % 
        % Preallocating space for each component of TF and for all stations
        ZXX = NaN(seg_eval,nts);           ZXY = NaN(seg_eval,nts);
        ZYX = NaN(seg_eval,nts);           ZYY = NaN(seg_eval,nts);
        TXZ = NaN(seg_eval,nts);           TYZ = NaN(seg_eval,nts);
        for a = 1:nts
            if act_chan{a}(1) && act_chan{a}(2)  
                % Extracting the best estimations for each component, for 
                % each station, for sorted Z
                ZXX(:,a) = Zxx(1 : seg_eval,a);
                ZXY(:,a) = Zxy(1 : seg_eval,a);
                ZYX(:,a) = Zyx(1 : seg_eval,a);
                ZYY(:,a) = Zyy(1 : seg_eval,a);
            end            
            if act_chan{a}(5)
                % Extracting the best estimations for each component, for 
                % each station, for sorted T
                TXZ(:,a) = Txz(1 : seg_eval,a);
                TYZ(:,a) = Tyz(1 : seg_eval,a);
            end
        end
        
        % Calculating weighted medians for transfer functions: Z, T and PT
        % for each time series (nts)
        for a = 1:nts
            % If impedances can be calculated (Ex and Ey are active)
            if act_chan{a}(1) && act_chan{a}(2)   
                % ------- Calculations for Impedances and Impedances Errors   
                [Zxx_wmed,Zxx_Err] = TF_weightederrors(ZXX(:,a),seg_eval,Wxx{a});
                [Zxy_wmed,Zxy_Err] = TF_weightederrors(ZXY(:,a),seg_eval,Wxy{a});
                [Zyx_wmed,Zyx_Err] = TF_weightederrors(ZYX(:,a),seg_eval,Wyx{a});
                [Zyy_wmed,Zyy_Err] = TF_weightederrors(ZYY(:,a),seg_eval,Wyy{a});
            % If impedances can not be calculated (Ex and Ey are inactive)
            % fill with NaN
            else
                Zxx_wmed = complex(0);      Zxx_Err = complex(0);
                Zxy_wmed = complex(0);      Zxy_Err = complex(0);
                Zyx_wmed = complex(0);      Zyx_Err = complex(0);
                Zyy_wmed = complex(0);      Zyy_Err = complex(0);                               
            end
            
            % If tipper can be calculated (Bz is/are active)
            if act_chan{a}(5)
                % ------------------- Calculations for Tipper and Tipper Errors                
                [txz_wmed,txz_Err] = TF_weightederrors(TXZ(:,a),seg_eval,Wxz{a});
                [tyz_wmed,tyz_Err] = TF_weightederrors(TYZ(:,a),seg_eval,Wyz{a});
            % If tipper can not be calculated (Bz is/are inactive) fill
            % with NaN
            else
                txz_wmed = complex(0);      txz_Err = complex(0);
                tyz_wmed = complex(0);      tyz_Err = complex(0);
            end
            
            % If impedances can be calculated, PT can also be calculated
            if ~isnan(Zxx_wmed)
                % ------------------------------- Calculations for Phase Tensor
                Z = [Zxx_wmed;Zyx_wmed;Zxy_wmed;Zyy_wmed];
                Z_err = [Zxx_Err;Zyx_Err;Zxy_Err;Zyy_Err];
                PT = FFproc_PT(Z);
                PT_err = FFproc_PT_err(Z,Z_err);
                PT.phimin_Err = PT_err(:,1);        PT.phimax_Err = PT_err(:,2);
                PT.alpha_Err = PT_err(:,4);         PT.beta_Err = PT_err(:,3);
                PT.theta_Err = PT_err(:,5);         PT.lambda_Err = PT_err(:,6);
            % If impedances can not be calculated, fill PT component with 
            % NaN
            else                
                PT.phi11 = complex(0);          PT.phi12 = complex(0);
                PT.phi21 = complex(0);          PT.phi22 = complex(0);
                PT.phimin = complex(0);         PT.phimax = complex(0);
                PT.alpha = complex(0);          PT.beta = complex(0);
                PT.theta = complex(0);          PT.lambda = complex(0);
                PT.phimin_Err = complex(0);     PT.phimax_Err = complex(0);
                PT.alpha_Err = complex(0);      PT.beta_Err = complex(0);
                PT.theta_Err = complex(0);      PT.lambda_Err = complex(0);
            end            
            
            % Saving transfer functions in results structure: Z 
            results(r).Zxx(a) =         Zxx_wmed;
            results(r).Zxy(a) =         Zxy_wmed;
            results(r).Zyx(a) =         Zyx_wmed;
            results(r).Zyy(a) =         Zyy_wmed;
            results(r).Zxx_Err(a) =     Zxx_Err;
            results(r).Zxy_Err(a) =     Zxy_Err;
            results(r).Zyx_Err(a) =     Zyx_Err;
            results(r).Zyy_Err(a) =     Zyy_Err;  
            % Saving transfer functions in results structure: T 
            results(r).txz(a) =         txz_wmed;
            results(r).tyz(a) =         tyz_wmed;
            results(r).txz_Err(a) =     txz_Err;
            results(r).tyz_Err(a) =     tyz_Err;
            % Saving transfer functions in results structure: PT
            results(r).phi11(a) =       PT.phi11;        
            results(r).phi12(a) =       PT.phi12;
            results(r).phi21(a) =       PT.phi21;        
            results(r).phi22(a) =       PT.phi22;
            results(r).phimax(a) =      PT.phimax;
            results(r).phimin(a) =      PT.phimin;      
            results(r).alpha(a) =       PT.alpha; 
            results(r).beta(a) =        PT.beta;
            results(r).theta(a) =       PT.theta;        
            results(r).lambda(a) =      PT.lambda;
            results(r).phimin_Err(a) =  PT.phimin_Err;
            results(r).phimax_Err(a) =  PT.phimax_Err;
            results(r).alpha_Err(a) =   PT.alpha_Err;
            results(r).beta_Err(a) =    PT.beta_Err;
            results(r).theta_Err(a) =   PT.theta_Err;
            results(r).lambda_Err(a) =  PT.lambda_Err;
        end
        
        % If number of time series (nts) is greater than 2, Horizontal
        % Magnetic Transfer Functions (HMTF) can be calculated.
        if nts > 1
            count_ch = 1;
            for i = 1:nts
                for j = 1:nts
                    if magref(i,1) == i && magref(i,2) == i && ...
                           magref(j,1) == j && magref(j,2) == j && ...
                           ~isequal(magref(i,:),magref(j,:))
                        count_ch = count_ch + 1;
                    end
                end
            end
        end
    end   

    % ------------------------------------------------------------------- %
    %                        Saving EigVal structure                      %
    % ------------------------------------------------------------------- % 
    % Eigevalue Index is now sorted and the indices are obtained
    [EI_sort,EI_ind] = sort(EI,'descend');

    % Populating EigVal structure with EI
    EigVal.EI0 =        EI0;
    EigVal.EI1 =        EI1;
    EigVal.EI2 =        EI2;
    EigVal.EI3 =        EI3;
    EigVal.EI4 =        EI4;
    EigVal.EI =         EI;
    EigVal.EI_sort =    EI_sort;       % Sorted Eigenvualue Index
    EigVal.EI_ind =     EI_ind;        % Sorted Indices
    EigVal.evpow =      evpow;

    % Populating EigVal structure with PCC
    EigVal.Zxx_pcorr =  cat(1,TF.rZxx);
    EigVal.Zxy_pcorr =  cat(1,TF.rZxy);
    EigVal.Zyx_pcorr =  cat(1,TF.rZyx);
    EigVal.Zyy_pcorr =  cat(1,TF.rZyy);
    EigVal.txz_pcorr =  cat(1,TF.rtxz);
    EigVal.tyz_pcorr =  cat(1,TF.rtyz);
    EigVal.rpow = rpow;

    % Populating EigVal structure with weighted EI * PCC
    EigVal.EIPC_Zxx = EIPC_Zxx;      EigVal.EIPC_Zxx_ind = EIPC_Zxx_ind;
    EigVal.EIPC_Zxy = EIPC_Zxy;      EigVal.EIPC_Zxy_ind = EIPC_Zxy_ind;
    EigVal.EIPC_Zyx = EIPC_Zyx;      EigVal.EIPC_Zyx_ind = EIPC_Zyx_ind;
    EigVal.EIPC_Zyy = EIPC_Zyy;      EigVal.EIPC_Zyy_ind = EIPC_Zyy_ind;
    EigVal.EIPC_txz = EIPC_txz;      EigVal.EIPC_txz_ind = EIPC_txz_ind;
    EigVal.EIPC_tyz = EIPC_tyz;      EigVal.EIPC_tyz_ind = EIPC_tyz_ind;    

    % Saving information for each segment
    EigVal.segments = segments;    
end

% results is a structure variable containing the following fields:
% Zxx, Zxy, Zyx, Zyy, Zxx_Err, Zxy_Err, Zyx_Err, Zyy_Err,
% txz, tyz, txz_Err, tyz_Err
% phi11, phi12, phi21, phi22, phimax, phimin, alpha, beta, theta,
% lambda, and errors...
% txx, txy, tyx, tyy, txx_Err, txy_Err, tyx_Err, tyy_Err, 
% ***_EVal (sorted eigenvalues for each component).