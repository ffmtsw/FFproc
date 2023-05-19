function EI = calcevindex(EV,ind,tfout,remmethod)

    if nargin < 3
        tfout = false;
    elseif nargin == 3
        remmethod = 'gesd';
    end

    % This is a test, we are trying to identify which weight is better: the
    % EV magnitude or the log10 of the EV magnitude (cc / 16mar2023)
    tflog = false;

    if all(isnan(EV))
        EV1 = ones(size(EV));
        EV2 = ones(size(EV));
        EV3 = ones(size(EV));
    else
        % Eigenvalue 1, 2 and 3
        if tflog
            EV1 = log10(abs(EV(:,1))); %#ok<*UNRCH> 
            EV2 = log10(abs(EV(:,2)));
            EV3 = log10(abs(EV(:,3)));
        else
            EV1 = abs(EV(:,1));
            EV2 = abs(EV(:,2));
            EV3 = abs(EV(:,3));
        end
    end

    % Removing possibles NaN or Inf and replacing by the smallest value
    % different from NaN or Inf
    EV1(isnan(EV1) | isinf(EV1)) = min(EV1(~isnan(EV1)|~isinf(EV1)));
    EV2(isnan(EV2) | isinf(EV2)) = min(EV2(~isnan(EV2)|~isinf(EV2)));
    EV3(isnan(EV3) | isinf(EV3)) = min(EV3(~isnan(EV3)|~isinf(EV3)));

    % Threshole for EV1 
    tfEV1 = EV1 < 1;
    EV1(tfEV1) = min(EV1(~isnan(EV1)|~isinf(EV1)));
    EV2(tfEV1) = min(EV1(~isnan(EV2)|~isinf(EV2)));
    EV3(tfEV1) = min(EV1(~isnan(EV3)|~isinf(EV3)));
    % Threshole for EV2 
    tfEV2 = EV2 < 1;
    EV1(tfEV2) = min(EV1(~isnan(EV1)|~isinf(EV1)));
    EV2(tfEV2) = min(EV1(~isnan(EV2)|~isinf(EV2)));
    EV3(tfEV2) = min(EV1(~isnan(EV3)|~isinf(EV3)));

    % Reduce weights to [0,1]?
    reduce = true;
    
    switch ind
        case 0  % Eigenvalue 0: |λ2|/|λ1|            
            EI0 = EV2./EV1;
            if tfout
                tf = isoutlier(EI0,remmethod);
                EI0(tf) = min(EI0(~tf));
            end
            if numel(EI0) > 1
                EI = EI0;
            else
                EI = 1;
            end
        case 1  % Eigenvalue 1: 1 - |λ3|/|λ2|            
            EI1 = 1-(EV3./EV2);
            if tfout
                tf = isoutlier(EI1,remmethod);
                EI1(tf) = min(EI1(~tf));
            end
            if numel(EI1) > 1
                EI = EI1;
            else
                EI = 1;
            end
        case 2  % Eigenvalue 2: 1/|λ3|            
            EI2 = 1./EV3;
            if tfout
                tf = isoutlier(EI2,remmethod);
                EI2(tf) = min(EI2(~tf));
            end
            if numel(EI2) > 1
                EI = EI2;
            else
                EI = 1;
            end
        case 3  % Eigenvalue 3: |λ2|/|λ3| - |λ3|^2   
            EI3 = (EV2./EV3) - EV3.^2;
            if tfout                
                tf = isoutlier(EI3,remmethod);
                EI3(tf) = min(EI3(~tf));
            end     
            if numel(EI3) > 1
                EI = EI3;
            else
                EI = 1;
            end
        case 4  % Eigenvalue 4: |λ2|/|λ3| - |λ1|/|λ2|           
            EI4 = (EV2./EV3) - (EV1./EV2);
            if tfout
                tf = isoutlier(EI4,remmethod);
                EI4(tf) = min(EI4(~tf));
            end    
            if numel(EI4) > 1
                EI = EI4;
            else
                EI = 1;
            end
    end

    if numel(EI) > 1
        if reduce
            EI = EI - min(EI);
            EI(EI < 0) = 0;
            EI = EI./max(EI);
        end
    end

    % Removing possibles NaN or Inf and replacing by the smallest value
    % different from NaN or Inf
    EI(isnan(EI) | isinf(EI)) = min(EI(~isnan(EI)|~isinf(EI)));

end