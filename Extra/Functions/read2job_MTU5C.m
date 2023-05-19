function job = read2job_MTU5C(root)

    if nargin == 0
        root = uigetdir(pwd,'Select folder containing MTU5-C measurements');
        cd(root)
    elseif isempty(root)
        root = uigetdir(pwd,'Select folder containing MTU5-C measurements');
        cd(root)
    else
        cd(root)
    end
    
    info = dir(root);
    info(strcmp(extractfield(info,'name'),'.') | strcmp(extractfield(info,'name'),'..')) = [];

    ts = [];
    for i = 1:numel(info)    
        try
            sitefolder = fullfile(root,info(i).name);
            cd(sitefolder)
        catch
        end

        try
            ext = {'.td_96k','.td_24k','.td_2400','.td_150','.td_30'};
            tf = false(numel(ext),1);
            for j = 1:numel(ext)
                tf(j) = ~isempty(dir(['**/*',ext{j}]));
            end
            % Removing not included sampling rates
            ext = ext(tf);
        catch
        end
        
        tsinput = [];
        for j = 1:numel(ext)
            try
                switch ext{j}
                    case '.td_96k'
                        sr = 96000; 
                    case '.td_24k'
                        sr = 24000;
                    case '.td_2400'
                        sr = 2400;
                    case '.td_150'
                        sr = 150;
                    case '.td_30'
                        sr = 30;
                end
                tsnew = read_MTU5C(sitefolder,sr,false,[],[]);
                tsinput = [tsinput,tsnew];
            catch
                tsinput = [];
            end
        end
    
        ts = [ts,tsinput];
    
    end
    
    %%
    for i = 1:numel(ts)
        job(i) = job_create;
        job(i) = jobconfig(job(i),ts(i),ts(i).device,ts(i).type);
    end

end