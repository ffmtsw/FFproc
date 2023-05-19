% function FFproc_K1_ANM is the 1/3 core function of FFproc_Process.m
% During this function the Noise Covariance Matrix (NCM), the Standar Noise
% Model (SNM), the Spectral Density Matrix (SDM), Eigenvalues, EigenVectors
% and Errors are calculated.
% version 1.0 / ?????2019 / ph
% version 2.0 / 12apr2020 / cc   Modified to be implemented within the new
%                                version of FFproc and FFMT (2.0).
%                                Name has been changed from
%                                processprocess_aut.m to EgbertHering.m
%                                No restriction on channels to be used.
%                                Many variables have been compressed in one
%                                to enhance the performance.
% version 2.1 / 17nov2020 / cc   Modified to improve performance by using
%                                arrays calculations instead of for
%                                loops. Active channels are included in
%                                order to perform Remote Reference
%                                calculations.
% version 2.2 / 24nov2020 / cc   Modified to remove diference between for
%                                and parfor (for serial or parallel 
%                                calculations).
% version 2.3 / 09feb2021 / cc   Option for saving Fourier Coefficients
%                                enabled. Information will be saved in FC
%                                structure.
% version 3.0 / 25aug2022 / cc   Bivariate linear regresssion included to
%                                estimate Transfer Functions. Option to
%                                detect outliers. Mahalanobis Distances has
%                                been removed.

% INPUT VARIABLES
% job:       global variable with information from FFproc_Parameters GUI
% Fmat:      matrix with Fourier coefficients (NCoefficientsxNTimeSegments)
% tchan:     number of total of channels (active and not active)
% act_chan:  logical vector with active channels to be used for each ts
% ntseg:     number of time segments
% parproc:   logical - to perform parallel or serial calculations

function [linreg,FC] = FFproc_K1_ANM(job,Fmat,act_chan,ntseg,parproc)
    
    warning('off','all')
    
    % Sampling rate
    sr = job.sr; %#ok<NASGU>
    
    % Number of time series per job
    nts = numel(job.ffts);
    
    % Total number of channels per job
    tchan = numel([act_chan{:}]);

    % Estimation approach
    if job.bivar
        approach = 'bivar';
    elseif job.ev
        approach = 'ev';
    end
    
    % Noise Model Estimator
    estimator = job.model;
         
    % Regression approach (Robust/Linear)
    reg = job.reg; 
    
    % Number of workers
    if parproc         
        ncores = job.ncores;
    else
        ncores = 0;
    end    

    % Outlier detection and removal
    tfout = job.outlier;
    remmethod = job.remmethod;
        
    % Output structure
    linreg = struct;
    
    % Indices for channels (ex, ey, bx, by, bz)
    ex = 1; ey = 2; bx = 3; by = 4; bz = 5;

    % Option for saving Fourier Coefficients
    FC_save = job.fc;

    % FC Structure
    FC = struct;

    % Magnetic reference channels
    magref = zeros(nts,2);
    for a = 1:nts
        magref(a,1) = job.ffts(a).rsbx;
        magref(a,2) = job.ffts(a).rsby;
    end
    
    % Extracting all Fourier coeficients matrices from all channels        
    Fdata = cat(1,Fmat);
    
    % Concatenating all Fourier coefficients in one matrix
    Fdata = cat(1,Fdata{:});    
    
    parfor (w = 1:ntseg,ncores)    
%     for w = 1:ntseg
        warning('off','all')

        % Each column is related to each time window
        F_seg = Fdata(:,w);
        % Reshaping Fourier coefficients to active channels size
        F_seg = reshape(F_seg,[],tchan);

        % --------------------------------------------------------------- %
        %                REMOVING SELECTED SPECTRAL LINES                 %
        % --------------------------------------------------------------- %
        % Removing Fourier coefficients == 0
        % coefficients with 0 come from unwanted frequencies
        F_seg(F_seg(:,ex) == 0 & F_seg(:,ey) == 0,:) = [];
        F_seg(F_seg(:,ex) == complex(0) & F_seg(:,ey) == complex(0),:) = [];
        F_seg(F_seg(:,bx) == 0 & F_seg(:,by) == 0,:) = [];
        F_seg(F_seg(:,bx) == complex(0) & F_seg(:,by) == complex(0),:) = [];

        % --------------------------------------------------------------- %
        %                          OUTLIER REMOVAL                        %
        % --------------------------------------------------------------- %
        if tfout
            % Number of Fourier Coefficients
            n_F = size(F_seg,1);
            switch remmethod                
                case 'maha'
                    try
                        % Test for outliers using Mahalanobis distances of the 
                        % observations using the robust estimates of the mean 
                        % and covariance. This is an iterative procedure with 3
                        % criteria to finish:
                        %  tfiter: no outliers detected
                        %  iter:   only 3 iterations
                        %  tf:     number outliers no less than 1/3 of the
                        %          number of input spectral lines
                        tfiter = true;
                        iter = 0;
                        tf = 0;
                        while tfiter && iter < 3 && tf < n_F*0.30
                            iter = iter + 1;
                            [~,~,~,tf_re] = robustcov(real(F_seg));
                            [~,~,~,tf_im] = robustcov(imag(F_seg));
                            tfmaha =  logical(tf_re + tf_im);
                            if sum(tfmaha) == 0
                                tfiter = false;
                            else
                                F_seg(tfmaha,:) = [];
                            end
                            tf = tf + sum(tfmaha);
                        end    
                    catch
                    end
                otherwise
                    % Test for outliers using MATLAB algorithms: mean,
                    % median, grubbs and gesd.
                    tf = false(size(F_seg));            
                    for o = 1:size(tf,2)
                        tf(:,o) = TF_downweight(F_seg(:,o),[],remmethod);
                    end
                    tf = logical(sum(tf,2));
                    if ~all(tf)
                        F_seg(tf,:) = [];
                    end
                    % Number of outliers
                    tf = sum(tf);
            end
        else
            tf = 0;
        end


        % ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| %
        % VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV %
        %                                                                 %
        %                       START OF CALCULATIONS                     %
        %                                                                 %
        % --------------------------------------------------------------- %
        % Input spectral coefficients
        F_in = F_seg;

        % Number of Fourier Coefficients
        n_F = size(F_in,1);

        % --------------------------------------------------------------- %
        %                   MULTIVARIATE LINEAR REGRESSION                %
        % --------------------------------------------------------------- %

        % Univariate or multivariate linear regression
        [fitted,NCM,rmse1,rmse2] = unimultivarlinreg(F_in,n_F,nts,act_chan,tchan,estimator,reg,approach);

        % ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| %
        % VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV %
        %                                                                 %
        %                        END OF CALCULATIONS                      %
        %                                                                 %
        % --------------------------------------------------------------- %

        % Saving regression results
        linreg(w).noutlier = tf;
        linreg(w).fitted = fitted;
        linreg(w).NCM = NCM;
        linreg(w).rmse1 = rmse1;
        linreg(w).rmse2 = rmse2;

        if FC_save
            FC(w).FC_Original = F_seg;
            FC(w).FC_Outlier = F_in;
            FC(w).FC_Fitted = fitted;
            FC(w).NCM = NCM;
            FC(w).RMS1 = rmse1;
            FC(w).RMS2 = rmse2;
        else
            FC(w).FC_Original = [];
            FC(w).FC_Outlier = [];
            FC(w).FC_Fitted = [];
            FC(w).NCM = [];
            FC(w).RMS1 = [];
            FC(w).RMS2 = [];            
        end

    end

end

% dataout is a structure variable containing the following fields:
% EVal, W1, W2, exbx_err, exby_err, eybx_err, eyby_err,
% bxbx_err and bybz_err.
% The number of elements of dataout are the number of time segments (numel(dataout) = ntseg)
% The dimensions of W1 and W2 are (Nactivechannels x Ntimeseries)
% The dimensions of errors are (Ntimeseries x 1)