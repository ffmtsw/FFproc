% function FFproc_datafromjob extracts information from Job Table
% this information can be from the original job or entered (typed) by the
% user. The new information is overwritten in job global variable.
% version 1.0 / 07abr2020 / cc

function FFproc_datafromjob(app,job)

    info = app.jobtable.Data;

    % Reasigning data to variables
    for i = 1:numel(job)
        % From stationtable
        job(i).sr = info{i,2}(1);
        job(i).ftarg = info{i,3}(1);
        job(i).maxf = info{i,4}(1);
        job(i).minf = info{i,5}(1);
        job(i).osc = info{i,6}(1);
        switch info{i,7}
            case 'Bivariate'
                job(i).bivar = true;
                job(i).ev = false;
            case 'EV Decomp'
                job(i).bivar = false;
                job(i).ev = true;
        end
        job(i).model = info{i,8}(1:end);
        job(i).reg = info{i,9}(1:end);
        job(i).outlier = info{i,10}(1);
        job(i).remmethod = info{i,11}(1:end);
        job(i).zrank = info{i,12}(1:end);
        job(i).trank = info{i,13}(1:end);
        job(i).ranking = info{i,14}(1);
        job(i).eind = info{i,15}(1:end);
        job(i).evpow = info{i,16}(1);
        job(i).rpow = info{i,17}(1);
        job(i).polyfit = info{i,18}(1);
        job(i).opt = info{i,19}(1);
        job(i).rnk = info{i,20}(1);
        job(i).raw = info{i,21}(1);
        job(i).fc = info{i,22}(1);
        job(i).ncores = info{i,23}(1);
    end

    app.job = job;
    
end