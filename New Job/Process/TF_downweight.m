function tf = TF_downweight(TF,nrem,outmet)

    % Number of Transfer Functions
    N = numel(TF);

    % If the ammount of data that is allowed to be removed is empty or not
    % parsed, as default, only up to 30% is allowed to be removed.
    if nargin == 1
        nrem = 0.3;
        % Number of maximum outliers to be detected
        outmax = ceil(N*nrem);
    else
        if isempty(nrem)
            nrem = 0.3;
        end
        % Number of maximum outliers to be detected
        outmax = ceil(N*nrem);
    end
    
    % Method to detect outliers
    if nargin == 3
        % Function contains detection method
    else
        outmet = 'gesd';
    end
    
    % Only if there are at least 2 spectral lines
    if N > 1

        switch outmet
            case 'maha'
                % Look for outliers based on the Mahalanobis distances
                iter = true;
                % Number of outliers
                niter = 0;
                while iter || niter < 2
                    % Ammount of outlieres from input data
                    N_out = find(isnan(TF));
                    if isempty(N_out)
                        N_out = 0;
                    end
                    % Only if the ammount of outliers is less than the upper limit
                    % allowed to remove outliers
                    if N_out/N < outmax
                        niter = niter + 1;
                        [~,~,~,tfr] = robustcov(real(TF));
                        [~,~,~,tfi] = robustcov(imag(TF));
                        % Outliers on real or imaginary part
                        out = logical(tfr + tfi);
                        % Asisgning NaN to TF outliers
                        TF(out) = complex(NaN,NaN);
                        if sum(tfr) == 0 & sum(tfi) == 0
                            iter = false;
                        else
                            niter = niter + 1;
                        end
                    else
                        iter = false;
                    end
                end
            case 'grubbs'
                tfr = isoutlier(real(TF),outmet);
                % Imaginary part
                tfi = isoutlier(imag(TF),outmet);    
                % Outliers on real or imaginary part
                out = logical(tfr + tfi);
                % Asisgning NaN to TF outliers
                TF(out) = complex(NaN,NaN);
            case 'gesd'
                % Real part
                tfr = isoutlier(real(TF),outmet,'MaxNumOutliers',outmax);
                % Imaginary part
                tfi = isoutlier(imag(TF),outmet,'MaxNumOutliers',outmax);    
                % Outliers on real or imaginary part
                out = logical(tfr + tfi);
                % Asisgning NaN to TF outliers
                TF(out) = complex(NaN,NaN);
            otherwise
                % Look for outliers based on the mean or median
                iter = true;
                % Number of outliers
                niter = 0;
                while iter    
                    % Ammount of outlieres from input data
                    N_out = find(isnan(TF));
                    if isempty(N_out)
                        N_out = 0;
                    end
        
                    % Only if the ammount of outliers is less than the upper limit
                    % allowed to remove outliers
                    if N_out/N < outmax
                        
                        % Real part
                        tfr = isoutlier(real(TF),outmet);
                        % Imaginary part
                        tfi = isoutlier(imag(TF),outmet);
            
                        % Outliers on real or imaginary part
                        out = logical(tfr + tfi);
        
                        % Asisgning NaN to TF outliers
                        TF(out) = complex(NaN,NaN);
            
                        if sum(tfr) == 0 & sum(tfi) == 0
                            iter = false;
                        else
                            niter = niter + 1;
                        end
                    else
                        iter = false;
                    end
                end
        end
    end

    % Indices of outliers
    tf = isnan(TF);

end