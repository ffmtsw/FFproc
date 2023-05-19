function [comp,error,select,fitcurve] = fitpoly2(per,arg,err,polyorder)

    %%%%%%%%%%%%%%%%%%%%%% NAN AND INFINITE REMOVAL %%%%%%%%%%%%%%%%%%%%%%%
    % Finding infinites and nans (Columns)
    inf_col = nansum(isinf(arg),1)>0; %#ok<*NANSUM> 
    if all(inf_col)
        ok_col = 1:size(arg,1);
    else
        nan_col = all(isnan(arg),1); 
        col = logical(inf_col + nan_col);
        ok_col = find(~col);
    end
    
    % Finding infinites and nans (Rows)
    for i = 1:size(arg,1)
        inf_row = isinf(arg(i,:));
        arg(i,inf_row) = median(arg(i,~inf_row));
        nan_row = isnan(arg(i,:));
        arg(i,nan_row) = median(arg(i,~nan_row));
    end
    inf_row = nansum(isinf(arg),2)>0;
    nan_row = all(isnan(arg),2);
    row = logical(inf_row + nan_row);
    ok_row = find(~row);
    
    % Selecting only data without inf and nan
    PER = per(ok_row);
    ARG = arg(ok_row,ok_col);    
    ERR = err(ok_row,ok_col);

    % Weighting Errors to be considered during robust fit
    W = ERR - min(ERR);
    W = 1 - W./max(W);
    W(isnan(W)) = 1;

    % Create matrix of repeated period vector
    PERMAT = repmat(PER,1,size(ARG,2));   

    %%%%%%%%%%%%%%%%%%%%%% PERIOD RANGE EXTENSION  %%%%%%%%%%%%%%%%%%%%%%%%
    tfextent = false;
    if tfextent
        % Number of extra decades
        bound = 1; %#ok<*UNRCH> 
        % Avoid edges problems by extending frequency range
        permin = min(PER);
        permax = max(PER);    
        perminplus = linspace(floor(permin) - bound,permin,10*bound);
        permaxplus = linspace(permax,ceil(permax) + bound,10*bound);
        PEREXT = [perminplus(1:end-1)';PER;permaxplus(2:end)'];        
    
        % Indices of limits: original arguments
        permin_id = find(PEREXT == permin);
        permax_id = find(PEREXT == permax);        
    
        % Argument extrapolation to new frequency range
        ARGEXT = interp1(PER,ARG,PEREXT,'linear','extrap'); 

        % Create matrix of repeated period vector
        PERMAT = repmat(PEREXT,1,size(ARGEXT,2));

        % Extending errors
        ERREXT = zeros(size(ARGEXT));
        ERREXT(permin_id:permax_id,:) = ERR;
        
        PER = PEREXT;
        ARG = ARGEXT;
        ERR = ERREXT;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% POLYNOMIAL FIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preallocate variables SSE and RSQ
    maxpoly = 9;
    RSQ = NaN(maxpoly,1);
    SSE = NaN(maxpoly,1);
    % Use all periods or extrema data
    tfextrema = false;
    % To start the while loop
    iter = false; 
    % Starting the counter
    n = 0;
    while ~iter        
        % Updating counter
        n = n + 1;        
        % GOF: Goodnes-of-Fitting
        if tfextrema
            PERMAT = PERMAT(:,1:2);
            ARGUMENT = [max(ARG,[],2);min(ARG,[],2)];
        else
            ARGUMENT = ARG;
        end
        [fitcurve,gof] = fit(PERMAT(:),ARGUMENT(:),['poly',num2str(n)],'Normalize','on','Robust','Bisquare','Weights',W(:));  
        % Error Sum of Squares (SSE)
        SSE(n) = gof.sse;
        % R-squared (R2)
        RSQ(n) = gof.rsquare;
        % Confidence Intervals
        cint = confint(fitcurve);
        cint_prod = cint(1,:).*cint(2,:);
       
        if any(cint_prod < 0) && RSQ(n) > 0.95 || n >= polyorder            
            % Last fitting over polynom order selected
            [~,id_polyord] = min(~isnan(SSE));
            [fitcurve,~] = fit(PERMAT(:),ARGUMENT(:),['poly',num2str(id_polyord)],'Normalize','on','Robust','Bisquare','Weights',W(:));
            iter = true;            
        end
    end    
    
    % Select id of the closest estimation to the fitted model
    SELECT = NaN(numel(PER),1);
    for f = 1:numel(PER)
        [~,SELECT(f)] = min(abs(ARG(f,:) - fitcurve(PER(f))));
    end   
    
    % Selecting the argument and error related to the closest estimation to 
    % the fitted model
    COMP = NaN(numel(PER),1);
    ERROR = NaN(numel(PER),1);
    for f = 1:numel(PER)        
        COMP(f) = ARG(f,SELECT(f));
        ERROR(f) = ERR(f,SELECT(f));
    end

    %%%%%%%%%%%%%%%%%%%%%%%% EXTRA DATA REMOVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Selection only within original frequency range
    if tfextent        
        COMP = COMP(permin_id:permax_id);
        ERROR = ERROR(permin_id:permax_id);
        SELECT = SELECT(permin_id:permax_id);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%% SAVING RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Saving output information
    comp = NaN(1,size(arg,1));    comp(ok_row) = COMP;
    error = NaN(1,size(arg,1));   error(ok_row) = ERROR;
    select = NaN(1,size(arg,1));  select(ok_row) = SELECT;    
    
end