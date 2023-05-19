function job = read2job_MTU5A(root)

    if nargin == 0
        root = uigetdir(pwd,'Select folder containing MTU5-A measurements');
        cd(root)
    elseif isempty(root)
        root = uigetdir(pwd,'Select folder containing MTU5-A measurements');
        cd(root)
    end
    
    info = dir(root);
    
    ts = [];
    for i = 1:numel(info)
    
        try
            cd(fullfile(root))
        catch
        end
        siteinfo = dir(fullfile(pwd,'*.TS*'));
        try
            tsinput = read_TS(siteinfo(i).folder,siteinfo(i).name,false,false);
        catch
            tsinput = [];
        end
    
        ts = [ts,tsinput];
    
    end
    
    %%
    for i = 1:numel(ts)
        job(i) = job_create;
        job(i) = jobconfig(job(i),ts(i),ts(i).device,ts(i).type);
    end

    if any(cat(1,job.sr) == 15)
        n = sum(cat(1,job.sr) == 15);        
        if n == 1
            id = find(cat(1,job.sr) == 15);
            nid = numel(job) + 1;
            job(nid) = job(id);

            % First frequencies
            job(id).osc = 100;
            job(id).maxf = 6;
            job(id).minf = 0.04;   

            % Last frequencies
            job(nid).osc = 25;
            job(nid).maxf = 0.06;
            job(nid).minf = 0.0004; 

        end
    end

end