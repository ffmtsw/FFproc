function [vecfit,residual,rmse] = modelregression(vec,mat,reg)
    
    switch reg
        case 'LSQ'
            model = mat\vec;
            vecfit = mat*model;
            % Residuals between predicted and observed data
            residual = vec - vecfit;
            rmse = rms(residual);
        case 'Robust'
            try
                [model,stats] = robustfit(mat,vec,'huber');
                vecfit = mat*model(2:end) + model(1);
                % Residuals between predicted and observed data
                residual = stats.resid;
                rmse = rms(stats.rstud);      
            catch
                % Due to few oints to perform robust estimation error could
                % appear. Try again with LSQ regression
                [vecfit,residual,rmse] = modelregression(vec,mat,'LSQ');
            end
        case 'LinRobust'
            % This method is not supported anymore, as fitlm is not working
            % for complex numbers since MATLAB R2020a release
            model = fitlm(mat,vec,'linear','RobustOpts','huber');
            vecfit = model.Fitted;
            % Residuals between predicted and observed data
            residual = model.Residuals.Raw;
            rmse = rms(model.Residuals.Studentized);
    end

end