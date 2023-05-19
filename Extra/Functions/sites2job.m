function sites2job

    % Load site 1
    [file1,path] = uigetfile('*.mat','Select Site 1');
    if file1 == 0
        return
    else
        load(fullfile(path,file1))
        job1 = job;
    end
    % Load site 2
    [file2,path] = uigetfile('*.mat','Select Site 2');
    if file2 == 0
        return
    else
        load(fullfile(path,file2))
        job2 = job;
    end
    
    % Sorting job 1 (if needed) based on sampling rate 1 
    [sr1,id1] = sort(cat(1,job1.sr),'descend');
    job1 = job1(id1);
    % Sorting job 2 (if needed) based on sampling rate 2 
    [sr2,id2] = sort(cat(1,job2.sr),'descend');
    job2 = job2(id2);
    
    
    % Detecting if sampling rates for each job matches
    sr_unique = sort(unique([sr1;sr2]),'descend');
    if all(ismember(sr1,sr2),2)
        tf = all([ismember(sr_unique,sr1),ismember(sr_unique,sr2)],2);
        sr_inc = sr_unique(tf);
    else        
        tf = all([ismember(sr_unique,sr1),ismember(sr_unique,sr2)],2);
        if sum(tf) == 0
            disp('Number of jobs in site 1 does not match with jobs in site 2')    
            return
        else
            sr_inc = sr_unique(tf);
        end
    end
    
    clear job
    % Pre-allocating job fields
    job(sum(tf)) = job_create;
    % Retrieving field names
    names = fieldnames(job);
    % Job information to be retrieved from job 1
    names = names(3:end);

    % Extracting time series from each job
    for j = 1:sum(tf)
    
            id1 = find(cat(1,job1.sr) == sr_inc(j));
            id2 = find(cat(1,job2.sr) == sr_inc(j));

            jobnew = job_create;
            % Site 1 name
            jobnew.name = job1(id1).name;          
            % Site 2 name
            jobnew.name(2) = job2(id2).name;


            for n = 1:numel(id1)
                % Info for site 1
                job1(id1(n)).ffts.ID = 1;
                job1(id1(n)).ffts.rsbx = 1;
                job1(id1(n)).ffts.rsby = 1;
                jobnew.ffts = job1(id1(n)).ffts;
            end
            
            for n = 1:numel(id1)
                % Info for site 2
                job2(id2(n)).ffts.ID = 2;
                job2(id2(n)).ffts.rsbx = 2;
                job2(id2(n)).ffts.rsby = 2;
                jobnew.ffts(2) = job2(id2(n)).ffts;
            end

            % Copying information from site 1 in new job
            for n = 1:numel(names)
                jobnew.(names{n}) = job1(id1(1)).(names{n});
            end
    
            job(j) = jobnew;
    
    end

    % Output file name: Site1-Site2
    filename = [strtok(file1,'.'),'-',strtok(file2,'.')];
    
    % Saving output file
    [file,path] = uiputfile([filename,'.mat'],'Save Job file:');
    if path ~= 0
        cd(path)
        save(fullfile(path,file),'job','-mat')
        try
            winopen(path)
        catch
        end
    else
        return
    end
    
end