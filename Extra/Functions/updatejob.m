function updatejob(app,jobold)

    if nargin == 1
        [file,path] = uigetfile('.mat','Select Job file you want to update');
        cd(path)
        jobin = load(fullfile(path,file));
        jobold = jobin.job;
    end

    oldnames = fieldnames(jobold);

    job(numel(jobold)) = job_create;
    newnames = fieldnames(job);
    tf = ismember(newnames,oldnames);

    % Channels of time series - information to be transferred old -> new
    trchannel = {'ex','ey','exchan','eychan','bxchan','bychan','bzchan','rsbx','rsby','actchan','exswap','eyswap','bxswap','byswap','calfact','RAP'};
  
    if nargin > 0
        dlg = uiprogressdlg(app.UIFigure,'Title','FFproc - Job Editor');
        gui = true;
    else
        gui = false;
    end

    n = 0; % Initial counter for failed conversions
    for j = 1:numel(jobold)
        if gui
            dlg.Message = ['Updating job ',num2str(j),' of ',num2str(numel(jobold))];
            dlg.Value = j/numel(jobold);
        end
        % Updating time series input file
        for a = 1:numel(jobold(j).ffts)
            try
                ts(a) = ts_input(jobold(j).ffts(a).path,jobold(j).ffts(a).file,false,[]);
                for c = 1:numel(trchannel)
                    try
                        ts(a).(trchannel{c}) = jobold(j).ffts(a).(trchannel{c});
                    catch ME
                        switch ME.message                        
                            case ['Unrecognized field name "',trchannel{c},'".']
                                switch trchannel{c}
                                    case 'bxchan'
                                        ts(a).(trchannel{c}) = true;
                                    case 'bychan'
                                        ts(a).(trchannel{c}) = true;
                                    otherwise
                                        disp(['Unrecognized field name "',trchannel{c},'".'])
                                end
                        end
                    end
                end
            catch ME
                txt{n + 1,1} = ['Job: ',num2str(j),' | Time Series: ',num2str(a),' | File was not found'];
                n = n + 1;
            end
            ts(a).ID = a;
        end
        job(j).ffts = ts;
        clear ts        
    end    
    if gui
        close(dlg)
    end

    if n > 0
        m = [{[num2str(n),' Job entries have errors']};...
            {'FFproc will not be able to process this jobs:'};...
            {''}];
        for i = 1:numel(n)
            m{end+1} = txt{i};
        end
        answer = uiconfirm(app.UIFigure,m,...
                                         'FFproc - Job Editor', ...
                                         'Options',{'Execute','Return'},'DefaultOption',2,'CancelOption',2,'Icon','Error');
        switch answer
            case 'Return'
                return
        end
    end



    for j = 1:numel(jobold)
        for n = 1:numel(tf)
            if tf(n)
                switch (newnames{n})
                    case 'ffts'
                    otherwise
                        job(j).(newnames{n}) = jobold(j).(newnames{n});
                end
            else
                switch newnames{n}
                    case 'osc'
                        job(j).osc = jobold(j).osz;
                    case 'eind'
                        job(j).eind = jobold(j).eig;
                    case 'ranking'
                        job(j).ranking = jobold(j).eigperc;
                    case 'outlier'
                        job(j).outlier = true;
                    case 'thres'
                        job(j).thres.active = jobold(j).thact;
                        job(j).thres.thres1 = jobold(j).thres1;
                        job(j).thres.thres2 = jobold(j).thres2;
                        job(j).thres.thres3 = jobold(j).thres3;
                end
            end
        end
    
        % New value added 01dic2022
        if job(j).outlier
            job(j).remmethod = 'gesd';
        else
            job(j).remmethod = 'N/A';
        end
        % New values added 15aug2022
        job(j).exclude = true;
        job(j).bivar = false;
        job(j).alpha = [];
        job(j).ev = true;
        job(j).eind = 'EV3';
        job(j).evpow = 1;
        job(j).zrank = 'EVPCC';
        job(j).trank = 'EVPCC';
        job(j).rpow = 1;
        job(j).polyfit = true; 
        job(j).rnk = true;
    end

    try
        cd(path)
        [file,path] = uiputfile('.mat','Save updated Job file',[strtok(file,'.'),'_updt.mat']);
    catch
        cd(pwd)
        [file,path] = uiputfile('.mat','Save updated Job file','_updt.mat');
    end
    
    cd(path)
    save(fullfile(path,file),'job')
    
end