function job = read2job_ATS(root)

    if nargin == 0
        root = uigetdir(pwd,'Select folder containing ATS measurements');
        cd(root)
    elseif isempty(root)
        root = uigetdir(pwd,'Select folder containing ATS measurements');
        cd(root)
    end
    
    info = dir(root);
    
    ts = [];
    for i = 1:numel(info)
    
        try
            cd(fullfile(root,info(i).name))
        catch
        end
        siteinfo = dir(fullfile(pwd,'*.ats'));
        try
            if isempty(ts)
                tsinput = read_ATS(siteinfo(1).folder,false,false,[]);
            else
                tsinput = read_ATS(siteinfo(1).folder,false,false,ts(1));
            end
        catch
            tsinput = [];
        end
    
        ts = [ts,tsinput];
    
    end
    
    %%
    for i = 1:numel(ts)
        job(i) = job_create;
        job(i) = jobconfig(job(i),ts(i),ts(i).device);
    end

end