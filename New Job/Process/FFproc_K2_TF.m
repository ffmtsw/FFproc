% function FFproc_K2_TF is  the 2/3 core functions of FFproc_Process.m
% This function is used to perform a first estimation of Impedances,
% Tipper components and Horizontal Magnetic Tranfer Functions.
% version 1.0 / ?????2019 / ph
% version 1.1 / ?????2019 / at
% version 2.0 / 12apr2020 / cc   Modified to be implemented within the new
%                                version of FFMT (2.0). Name has changed
%                                from TF_aut.
% version 2.0 / 12apr2020 / cc   Modified to avoid tact_chan and those     
%                                variables that were used to count
%                                channels. Instead of that, cell arrays are
%                                used to improve the performance and the
%                                counter_ch and counter_ts variables are
%                                used to assign correct channels. HMTF and
%                                RR processing were added.
% version 2.2 / 24nov2020 / cc   Modified to remove diference between for
%                                and parfor (for serial or parallel 
%                                calculations).
% version 2.3 / 25nov2020 / cc   Check for coil polarity was added (bxswap
%                                and byswap)

% NCM: Covariance Matrix of the noise (Sigma_N)

function [TF,FC] = FFproc_K2_TF(job,linreg,act_chan,ntseg,FC,parproc)

    warning('off','all')
           
    % Number of time series per job
    nts = numel(job.ffts);
    
    % Number of workers
    if parproc         
        ncores = job.ncores;
    else
        ncores = 0;
    end 

    % Bivariate estimation of Z and T
    bivar = job.bivar;
    
    % Option for saving Fourier Coefficients
    FC_save = job.fc;
    
    % Magnetic reference channels
    magref = zeros(numel(job.ffts),2);
    for a = 1:nts
        magref(a,1) = job.ffts(a).rsbx;
        magref(a,2) = job.ffts(a).rsby;
    end
    
    % Indices for channels (ex, ey, bx, by, bz)
    ex = 1; ey = 2; bx = 3; by = 4; bz = 5;
    
    TF = struct;

    parfor (w = 1:ntseg,ncores) 
%     for w = 1:ntseg

        % Number of outliers
        TF(w).noutlier = linreg(w).noutlier;

        % Clean (fitted) Fourier Coefficients
        F_out = linreg(w).fitted;

        % Number of Fourier Coefficients
        n_F = size(F_out,1);

        % Noise Covariance Matrix: calculated from the residuals of the
        % multivariate linear regression between all cleaned data channels
        NCM = linreg(w).NCM;        

        if bivar

            % ----------------------------------------------------------- %
            %                                                             %
            %                     BIVARIATE ESTIMATION                    %
            %                                                             %
            % ----------------------------------------------------------- %
            % Bivariate estimation of Impedance Tensor (Z),Tipper vector 
            % (T) and partial coherences rZ and rT
            
            % Number of channels
            nchan = size(cell2mat(act_chan),2);
            for a = 1:nts
                chjump = nchan*(a - 1);
                rchjumpx = nchan*(magref(a,1)-1);
                rchjumpy = nchan*(magref(a,2)-1);
                if act_chan{a}(1) && act_chan{a}(2) && ~act_chan{a}(3) && ~act_chan{a}(4)
                    F_out(:,bx + chjump) = F_out(:,bx + rchjumpx);
                    F_out(:,by + chjump) = F_out(:,by + rchjumpy);
                end
            end

            [Zxx,Zxy,Zyx,Zyy,txz,tyz] = bivariate(n_F,nts,F_out,act_chan,magref);

            % Partial Coherences rZ and rT
            try
                [rZxx,rZxy,rZyx,rZyy,rtxz,rtyz] = pcoherences(n_F,nts,F_out,act_chan,magref);
            catch
                rZxx = 0;   rZxy = 0;
                rZyx = 0;   rZyy = 0; 
                rtxz = 0;   rtyz = 0;
            end

            % Saving Eigenvalues
            TF(w).EVal = complex(NaN,NaN);

            for a = 1:nts
                % Impedance Tensor
                TF(w).Zxx(a) = -Zxx(a);
                TF(w).Zxy(a) = -Zxy(a);
                TF(w).Zyx(a) = -Zyx(a);
                TF(w).Zyy(a) = -Zyy(a);
                % Tipper vector
                TF(w).txz(a) = txz(a);
                TF(w).tyz(a) = tyz(a);
                % Partial Coherencies
                TF(w).rZxx(a) = rZxx(a);
                TF(w).rZxy(a) = rZxy(a);
                TF(w).rZyx(a) = rZyx(a);
                TF(w).rZyy(a) = rZyy(a);
                TF(w).rtxz(a) = rtxz(a);
                TF(w).rtyz(a) = rtyz(a);
            end          

            % Spectral density matrix is empty as cannot be calculated
            % following a bivariate approach
            SDM = [];
            SDMR = [];

        else
            % ----------------------------------------------------------- %
            %                                                             %
            %                    EIGENVALUE DECOMPOSITION                 %
            %                                                             %
            % ----------------------------------------------------------- %            

            % Spectral Density Matrix (SDM) 
            % -> Equals the signal covariance matrix (Hering, 2019)
            SDM = cov(F_out(:,[act_chan{:}]));
    
            % Robust-Scaled Spectral Density Matrix
            SDMR = NCM^(-1/2)* SDM * NCM^(-1/2);
    
            % Solving the eigenvalue problem and calculating eigenvectors 
            % and eigenvalues (u)
            % SDMR (u) = λ (u)
            if any(isnan(SDMR(:))) || any(isinf(SDMR(:)))
                I = complex(0,0);
                eigvec = repmat(I,size(SDMR));
                eigval = repmat(I,size(SDMR));
            else 
                [eigvec,eigval] = eig(SDMR);
            end
            EVal = diag(eigval);
            [~,vec_ind] = sort(diag(eigval),'descend');
            EVec = eigvec(:,[vec_ind(1),vec_ind(2)]);
            
            % Calculation of u (nchan*2 matrix, eigenvectors associated with 
            % the two largest eigenvalues from Evec (largest variances)
            u = NCM^(1/2) * EVec * (ctranspose(EVec) * NCM^(-1) * EVec) ^(-1);        
    
            % Re-asigning empty values as 0, if a channel is not active.                  
            allchan = [act_chan{:}];                                 
            
            % Re-asigning empty values as 0, if a channel is not active. 
            U = complex(zeros(numel(allchan),2));

            % U is k-dimensional complex matrix which ideally represent the
            % magnetic and electric fields that would be observed at all
            % sites for idealized quasi-uniform magnetic sources, linearly
            % polarized (Egbert, 1997).
            U(allchan,:) = u;
            
            % Saving sorted Eigenvalues (Largest to smallest)
            TF(w).EVal = EVal(vec_ind);

            % Arranging W1 in cell arrays with size (Nstations x 1) for
            % each estimated time segment            
            W = cell(1,nts);
            count_ch = 1;
            for a = 1:nts

                for b = 1:numel(act_chan{a}) %#ok<*PFBNS>
                    W{1,a}(b,:) = U(count_ch,:);
                    count_ch = count_ch + 1;
                end
            end

            for a = 1:nts                
                % Impedance Tensor estimation (Z) just if electric channels
                % are active. Remote Reference is used if the magnetic
                % reference channel index is different than the current time
                % series to be evaluated.
                if act_chan{a}(1) && act_chan{a}(2)                    
                    Z = [W{a}(ex,1),W{a}(ex,2);...
                         W{a}(ey,1),W{a}(ey,2)]*...
                        [W{magref(a,1)}(bx,1), W{magref(a,1)}(bx,2);...
                         W{magref(a,2)}(by,1), W{magref(a,2)}(by,2)]^(-1);
    
                    % Impedances components (complex conjugate)
                    Zxx = -(Z(1,1)');   Zxy = -(Z(1,2)');
                    Zyx = -(Z(2,1)');   Zyy = -(Z(2,2)');               
                    
                    % Saving Impedances in TF structure
                    TF(w).Zxx(a) = Zxx;
                    TF(w).Zxy(a) = Zxy;
                    TF(w).Zyx(a) = Zyx;
                    TF(w).Zyy(a) = Zyy;
                else
                    TF(w).Zxx(a) = complex(0);
                    TF(w).Zxy(a) = complex(0);
                    TF(w).Zyx(a) = complex(0);
                    TF(w).Zyy(a) = complex(0);
                end

                % Tipper components estimation just if vertical magnetic
                % field is active. Remote Reference is used if the magnetic
                % referenc channel index is different from the current time
                % series to be evaluated.
                if act_chan{a}(5)
                    T = [W{a}(bz,1),W{a}(bz,2)]*[W{magref(a,1)}(bx,1), W{magref(a,1)}(bx,2);...
                                                 W{magref(a,2)}(by,1), W{magref(a,2)}(by,2)]^(-1);                    
                    % Tipper components (complex conjugate)
                    TF(w).txz(a) = T(1)'; 
                    TF(w).tyz(a) = T(2)';
                else
                    TF(w).txz(a) = complex(0); 
                    TF(w).tyz(a) = complex(0);
                end             
            end         

             % Partial Coherences rZ and rT
            try
                [rZxx,rZxy,rZyx,rZyy,rtxz,rtyz] = pcoherences(n_F,nts,F_out,act_chan,magref);
            catch
                rZxx = 0;   rZxy = 0;
                rZyx = 0;   rZyy = 0; 
                rtxz = 0;   rtyz = 0;
            end

            for a = 1:nts
                % Partial Coherencies
                TF(w).rZxx(a) = rZxx(a);
                TF(w).rZxy(a) = rZxy(a);
                TF(w).rZyx(a) = rZyx(a);
                TF(w).rZyy(a) = rZyy(a);
                TF(w).rtxz(a) = rtxz(a);
                TF(w).rtyz(a) = rtyz(a);
            end

        end

        % Horizontal Magnetic Transfer Functions (HMTF) can be
        % estimated only if multivariate processing is used.
        if nts > 1            
            count_ch = 1;
            for a = 1:nts
                if ~bivar
                    for b = 1:nts
                        if magref(a,1) == a && magref(a,2) == a && ...
                           magref(b,1) == b && magref(b,2) == b && ...
                           ~isequal(magref(a,:),magref(b,:))
                            HMTF = [W{a}(bx,1), W{a}(bx,2);...
                                    W{a}(by,1), W{a}(by,2)]*...
                                   [W{b}(bx,1), W{b}(bx,2);...
                                    W{b}(by,1), W{b}(by,2)]^(-1);
                            % Saving HMTF in TF structure
                            TF(w).txx(1,count_ch) = HMTF(1,1);    TF(w).txy(1,count_ch) = HMTF(1,2);
                            TF(w).tyx(1,count_ch) = HMTF(2,1);    TF(w).tyy(1,count_ch) = HMTF(2,2); 
                            TF(w).t_id(:,count_ch) = [a,b];
                            count_ch = count_ch + 1;
                        end
                    end
                else
                TF(w).txx = [];     TF(w).txy = [];
                TF(w).tyx = [];     TF(w).tyy = [];
                TF(w).t_id = [];
                end 
            end
        else
            TF(w).txx = [];     TF(w).txy = [];
            TF(w).tyx = [];     TF(w).tyy = [];
            TF(w).t_id = [];
        end

        % Saving Fourier Coefficients and Regression output
        if FC_save
            FC(w).SDM = SDM;
            FC(w).SDMR = SDMR;
        else
            FC(w).SDM = [];
            FC(w).SDMR = [];
        end

    end       
        
end

% TF is a structure variable containing the following fields:
% Zxx, Zxy, Zyx, Zyy, txz, tyz, txz_Err, tyz_Err, txx, txy, tyx, tyy, t_id
% The number of elements of dataout are the number of time segments (numel(TF) = ntseg)
% The dimensions of each variable are (Ntimeseries x 1)