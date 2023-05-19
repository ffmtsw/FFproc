% Outlier-Removal of normally distributed data according to F-test.

% version 1.1 / 24aug2022 / aj
% version 1.2 / 24aug2022 / cc

% INPUT
% X:        Data (real (e.g. time series) or complex (e.g. Fourier Coefficients))
% W:        Weights (true = data is used, false = data is skipped)
% alpha:    Significance Level (alpha = [], the program calculates alpha)
% f_alpha:  Factor for increasing alpha (only if alpha is empty)

% OUTPUT
% tf:       logical vector (true = outlier detected)
% tfout:    index of detected outliers (from max to min)

function [tf,tfout] = remoutlier(X,alpha,f_alpha)

    cond = true;

    % Initial outlier counter
    o = 0;

    % Initial detected outliers (all NaN)
    tfout = NaN(size(X));

    % Initial index for outliers (all false)
    tf = false(size(X));

    % Initial state for input data (all true)
    W = true(size(X));
    
    % Absolute value of input data (Chi-squeared distribution)
    X2 = abs(X).^2;

    while cond

        % Number of active data elements
        N = sum(W);    

        % Check if the number is real or complex
        if isreal(X(1))
            fact = 1;   % real
        else
            fact = 2;   % complex
        end

        [X2_max,X2_max_ind] = max(X2.*W);
        if isempty(alpha) || alpha == 0
            alpha = 1 - (1/(f_alpha*N*fact));
        end

        % Inverse Cumulated Distribution Function
        S_alpha = icdf('F',alpha,fact,N*fact);

        Xa = X2_max/fact;                           % fact*1 degree of freedom
        Xb = (sum(X2.*W) - X2_max)/((N-1)*fact);    % fact*(N-1) degrees of freedom
        q = Xa/Xb;

        if q > S_alpha
            o = o + 1;
            tf(X2_max_ind) = true;          % Outlier detected
            W(X2_max_ind) = false;          % Data element deactivated
            tfout(o) = X2_max_ind;          % Number of outliers
        else
            cond = false;
        end

    end

    % Removing extra NaN
    tfout = tfout(~isnan(tfout));

end