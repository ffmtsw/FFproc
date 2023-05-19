% Function unimultivarlinreg performs the multivariate linear regression
% between:
% 
% - Channels of each station (single site)
% - Channels between two or more stations (multisite)
%
% For single-site, we try to reduce the impact of coherent noise in each
% channel, by performing the regression of the desired channel
% (observations) and the remaining channels (predictors). This would reduce
% the effect of possibles electrochemical reactions, different for each 
% electrode, and alterarions in single magnetic sensors.
% (man-made noise sources)
% For multi-site, we try to reduce the impact of coherent noise between
% a selected channel (observations) and the available channels from the 
% other stations (predictors). This would reduce the effect of antropic
% noise sources.
% The linear regression is performed twice.
% 
% Input variables:
%   data:       matrix with all channels
%   N:          number of spectral lines
%   nts:        number of MT stations
%   act_chan:   LOGICAL index for channels to be used during regression
%   tchan       number of all channels avaiable from all MT stations
%   estimator:  linear regression approach (Standar, Advanced)
%   regmethod:  linear regression method to be used: LSQ, Robust,
%               LinRobust. LinRobust is not working anymore since MATLAB
%               2020b release. However, we expect that will be reintroduced
%               soon and therefore, not removed as regression method.
%
% Output variables:
%   datafit:    predicted data by the model regression
%   res:        residuals for each channel
%   err1:       errors for first regression
%   err2:       errors for second regression

function [datafit,res,err1,err2] = unimultivarlinreg(data,N,nts,act_chan,tchan,estimator,regmethod,approach)

    err1 = zeros(1,tchan);
    err2 = zeros(1,tchan);

    if nts == 1
        estimator = 'Standard';
    end

    switch estimator
        case 'Standard'
            % ----------------------------------------------------------- %
            %                      FIRST REGRESSION                       %
            %                     (Signal prediction)                     %
            % ----------------------------------------------------------- %
            D = data;
            active = cat(2,act_chan{:});
            % Predicted data (empty space)
            predicted = zeros(N,numel(active));
            for b = 1:tchan % Number of channels
                if active(b)
                    % Observations Vector
                    observations = D(:,b);
                    % Predictors Matrix
                    predictors = D(:,1:numel(active)~=b & active);
                    % Predicted data by linear regression (observations vs
                    % predictors), no residuals are needed so far.
                    [predicted(:,b),~,err1(1,b)] = modelregression(observations,predictors,regmethod);
                end
            end
            datafit = predicted;

            % ----------------------------------------------------------- %
            %                      SECOND REGRESSION                      %
            %         (Signal prediction and residuals calculation)       %
            % ----------------------------------------------------------- %
            D = datafit;
            % Predicted data (empty space)
            predicted = zeros(N,numel(active));
            % Residuals Matrix
            residual = zeros(N,numel(active));
            for b = 1:tchan % Number of channels
                if active(b)
                    % Observations Vector
                    observations = D(:,b);
                    % Predictors Matrix
                    predictors = D(:,1:numel(active)~=b & active);
                    % Predicted data by linear regression (observations vs predictors)
                    [predicted(:,b),residual(:,b),err2(1,b)] = modelregression(observations,predictors,regmethod);
                end
            end            
            % Covariance of residuals
            res = diag(var(residual(:,active)));

            % Predicted data (cleaned channels)
            switch approach
                case 'bivar'
                    datafit = predicted;
            end

        case 'Advanced'
            % ----------------------------------------------------------- %
            %                      FIRST REGRESSION                       %
            %                     (Signal prediction)                     %
            %                       Data estimation                       %
            % ----------------------------------------------------------- %
            count_ts = 1;   % Counter for time series
            % Allocating cell for Z
            D = cell(nts,numel(cat(2,act_chan{:}))/nts);
            for a = 1:nts
                for b = 1:numel(act_chan{a}) 
                    if act_chan{a}(b)
                        D{a,b} = data(:,count_ts);
                    else
                        D{a,b} = {};
                    end
                    count_ts = count_ts + 1;
                end
            end

            count_ts = 1;       % Counter for time series
            count_ch = 1;       % Counter for active channels
            datafit = [];       % Matrix to save predicted models

            % Calculation of first regression (predicted signals)
            for a = nts:-1:1  
                % Predicted data (empty space)
                predicted = zeros(N,numel(act_chan{count_ts}));
                % Predictors Matrix
                predictors = cat(2,D{a,act_chan{a}});                    
                for b = 1:numel(act_chan{count_ts})
                    if act_chan{count_ts}(b)
                        % Observations Vector
                        observations = D{count_ts,b};
                        % Predicted data by linear regression (observations
                        % vs predictors), no residuals are needed so far.
                        [predicted(:,b),~,err1(1,count_ch)] = modelregression(observations,predictors,regmethod);
                        % Jumping to next active channel
                        count_ch = count_ch + 1; 
                    end                        
                end
                datafit = [datafit,predicted];
                % Jumping to next time series
                count_ts = count_ts + 1;
            end

            % ----------------------------------------------------------- %
            %                      SECOND REGRESSION                      %
            %         (Signal prediction and residuals calculation)       %
            %              Noise Covariance Matrix formulation            %
            % ----------------------------------------------------------- %
            count_ts = 1;       % Counter for time series
            for a = 1:nts
                for b = 1:numel(act_chan{a})
                    if act_chan{a}(b)
                        D{a,b} = datafit(:,count_ts);
                    else
                        D{a,b} = {};
                    end
                    count_ts = count_ts + 1;
                end
            end
            
            count_ts = 1;       % Counter for time series
            count_ch = 1;       % Counter for active channels
            data2_fit = [];     % Space for fitted data
            res = [];           % Space for residuals

            % Calculations for second regresion (predicted signals and residuals)
            for a = nts:-1:1
                % Predicted data (empty space)
                predicted = zeros(N,numel(act_chan{count_ts}));
                % Predictors Matrix
                predictors = cat(2,D{a,act_chan{a}});
                % Residuals Matrix
                residual = zeros(N,numel(act_chan{count_ts}));                                    
                for b = 1:numel(act_chan{count_ts})                                             
                    if act_chan{count_ts}(b)  
                        % Observations Vector   
                        observations = D{count_ts,b};
                        % Predicted data by linear regression (observations vs predictors)
                        [predicted(:,b),residual(:,b),err2(1,count_ch)] = modelregression(observations,predictors,regmethod); 
                        % Jumping to next active channel
                        count_ch = count_ch + 1; 
                    end                        
                end
                data2_fit = [data2_fit,predicted];
                % Covariance of residuals
                covres = cov(residual);
                % Saving residuals
                res(size(res,1)+1 : count_ch-1,size(res,2)+1 : count_ch-1) = covres(act_chan{count_ts},act_chan{count_ts});
                % Jumping to next time series
                count_ts = count_ts + 1; 
            end    

            % Predicted data (cleaned channels)
            switch approach
                case 'bivar'
                    datafit = data2_fit;
            end
    end

end