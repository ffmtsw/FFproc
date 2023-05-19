% function FFproc_stationtable extracts information from job global
% variable to be shown in the Job table.
% version 1.0 / 07abr2020 / cc
% version 1.1 / 14jun2020 / cc       Option to save RAW Files 

function FFproc_filljobtable(app,job)

    for i = 1:numel(job)

        nj{i,1} = i;
        nts{i,1} = numel(job(i).ffts);
        sr{i,1} = min(cat(1,job(i).ffts.sr));
        ftarg{i,1} = job(i).ftarg;
        maxf{i,1} = job(i).maxf;
        minf{i,1} = job(i).minf;
        ftarg{i,1} = job(i).ftarg;
        osc{i,1} = job(i).osc;

        approachopt = ({'Bivariate','EV Decomp'});
        if job(i).bivar
            approach{i,1} = 'Bivariate';
        elseif job(i).ev
            approach{i,1} = 'EV Decomp';
        end
        app.jobtable.ColumnFormat{1,7} = approachopt;
        
        mdopt = ({'Standard','Advanced'});
        model{i,1} = job(i).model;
        app.jobtable.ColumnFormat{1,8} = mdopt;

        % rgopt = ({'LSQ','Robust','PCA'});
        rgopt = ({'LSQ','Robust'});
        reg{i,1} = job(i).reg;
        app.jobtable.ColumnFormat{1,9} = rgopt;

        outlier{i,1} = job(i).outlier;

        remopt = ({'mean','median','grubbs','gesd','maha'});
        remmethod{i,1} = job(i).remmethod;
        app.jobtable.ColumnFormat{1,11} = remopt;

        zrankopt = ({'EV','PCC','EVPCC'});
        zrank{i,1} = job(i).zrank;
        app.jobtable.ColumnFormat{1,12} = zrankopt;

        trankopt = ({'EV','PCC','EVPCC'});
        trank{i,1} = job(i).trank;
        app.jobtable.ColumnFormat{1,13} = trankopt;

        ranking{i,1} = job(i).ranking;
    
        eind{i,1} = job(i).eind;

        evpow{i,1} = job(i).evpow;

        rpow{i,1} = job(i).rpow;
            
        polyfit{i,1} = job(i).polyfit;
        rnk{i,1} = job(i).rnk;
        opt{i,1} = job(i).opt;
        raw{i,1} = job(i).raw;
        fc{i,1} = job(i).fc;

        ncores{i,1} = job(i).ncores;
    end
    
    % Writing on table
    fill = [nts,sr,ftarg,maxf,minf,osc,...
            approach,model,reg,outlier,remmethod,...
            zrank,trank,ranking,eind,evpow,rpow,...
            polyfit,rnk,...
            opt,raw,fc,...
            ncores];

    app.jobtable.Data = fill;
    
end