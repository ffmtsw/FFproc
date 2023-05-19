function FFproc_mergeOPTandRAW(varargin)

    if nargin == 0 || isempty(varargin)
        path = uigetdir(pwd,'Select folder with processing results');
        if path == 0
            return
        end    
        tfgui = false;
    elseif nargin == 1
        path = varargin{1};
        tfgui = false;
    else
        app = varargin{1};
        path = varargin{2};
        tfgui = true;
        dlg = uiprogressdlg(app.FFproc,'Title','FFproc','ShowPercentage','on','Cancelable','on');
    end

    [~,foldername,~] = fileparts(path);
    cd(fullfile(path))
    infofolder = dir('*mat');
     
    if tfgui
        dlg.Message = {'Performing Auto-Combine in folder:';foldername};
    end
   
    % Split name
    def = struct;
    for j = 1:numel(infofolder)
        namesplit = strsplit(infofolder(j).name,'_');
        switch namesplit{1}
            case {'ANM','STD'}
                def(j).name = fullfile(infofolder(j).folder,infofolder(j).name);
                def(j).model = namesplit{1};
                def(j).reg = namesplit{2};
                %
                tfalpha = isstrprop(namesplit{3},'alpha');
                def(j).osc = str2double(namesplit{3}(~tfalpha));
                %
                freq = strsplit(namesplit{4},'-');
                freqmax = freq{1};
                freqmin = freq{2};

                tf = isstrprop(freqmax,'digit');
                freqmaxnum = str2double(freqmax(tf));
                tf = isstrprop(freqmax,'alpha');
                freqmaxunit = freqmax(tf);
                if strcmp(freqmaxunit,'s')
                    freqmaxnum = 1/freqmaxnum;
                end
                tf = isstrprop(freqmin,'digit');
                freqminnum = str2double(freqmin(tf));
                tf = isstrprop(freqmin,'alpha');
                freqminunit = freqmin(tf);
                if strcmp(freqminunit,'s')
                    freqminnum = 1/freqminnum;
                end
                def(j).freqlabel = namesplit{4};
                def(j).freqmax = freqmaxnum;
                def(j).freqmin = freqminnum;
                if numel(namesplit) == 10
                    % If RR was deployed
                    def(j).rr = false;
                    % Type of file, either OPT or RAW
                    type = strtok(namesplit{10},'.');
                    def(j).file = type(1:3);
                    def(j).perc = str2double(type(4:end-1));
                    def(j).filelabel = namesplit{10};
                elseif numel(namesplit) == 11
                    % If RR was deployed
                    def(j).rr = true;
                    % Type of file, either OPT or RAW
                    type = strtok(namesplit{11},'.');
                    def(j).file = type(1:3);
                    def(j).perc = str2double(type(4:end-1));
                    def(j).filelabel = namesplit{11};
                end
            otherwise
        end
    end

    % Let's go again, let's start wit loops
    typeopt = {'FIT','SLD','ALL'};    
    modelopt = {'STD','ANM'};
    regopt = {'LSQ','ROB','PCA'};
    
    for t = 1:numel(typeopt)
        % Find different type of files ('FIT','SLD','ALL')
        tf_type = all(ismember(cat(1,def.file),typeopt{t}),2);
        for m = 1:numel(modelopt)
            % Find different type of noise models ('STD','ANM')
            tf_model = all(ismember(cat(1,def.model),modelopt{m}),2);              
            for r = 1:numel(regopt)
                % Find different type of Linear regression approach ('LSQ','ROB','PCA')
                tf_reg = all(ismember(cat(1,def.reg),regopt{r}),2);
                % Logical index of files of the same type
                tf = all([tf_type,tf_model,tf_reg],2);
                if sum(tf) > 0
                    % Extracting info of selected sites
                    str = def(tf);
                    if ~isempty(str)
                        % Loading and concatenating each file
                        tmpFFsave = [];
                        tmpjob = [];
                        for f = 1:numel(str)                            
                            load(str(f).name,'FFsave','job');
                            if f == 1
                                N = numel(job.ffts);
                            end
                            tmpFFsave = [tmpFFsave,FFsave];
                            tmpjob = [tmpjob,job];
                        end        
                        % Merge for each site
                        site = cat(1,tmpFFsave.site);
                        uniquesite = unique(site);
                        for u = 1:numel(uniquesite)
                            ind = ismember(site,uniquesite(u));
                            ffsave = tmpFFsave(ind);
                            % Sorting sampling rates from EV structure
                            [~,id] = sort(cat(1,tmpjob.sr),'descend');
                            ffsave = ffsave(id);
                            try
                                mergeFFsave(ffsave(cat(1,str.rr)),uniquesite(u),modelopt(m),regopt(r),typeopt(t),true);
                            catch
                            end
                            try
                                mergeFFsave(ffsave(~cat(1,str.rr)),uniquesite(u),modelopt(m),regopt(r),typeopt(t),false);
                            catch
                            end
                        end                        
                    end
                else
                    continue
                end           

            end
        end
    end
    infofolder = dir('*mat');
    for u = 1:numel(uniquesite)
        MT = [];
        for i = 1:numel(infofolder)
            namesplit = strsplit(infofolder(i).name,'_');
            if matches(uniquesite{u},namesplit{1})
                load(fullfile(infofolder(i).folder,infofolder(i).name)) %#ok<*LOAD> 
                MT = [MT,mt]; %#ok<AGROW> 
                delete(fullfile(infofolder(i).folder,infofolder(i).name))
            end
        end
        mt = MT;
        filesave = ['AUTCMB_',uniquesite{u},'.mat'];
        if isfile(filesave)
            delete(filesave)
        end
        save(filesave,'mt')
        clear mt
    end
    if tfgui
        uialert(app.FFproc,'Auto-Combine process: Completed','FFproc','Icon','success')
    end
end

%% Function to create FFsave files
function mergeFFsave(ffsave,name,model,reg,type,rr)
    
    [freq,ind,~] = unique(cat(1,ffsave.freq),'stable');
    lonlat = ffsave(1).lonlat;
    [x,y,~] = deg2utm(lonlat(1),lonlat(2));
    z = ffsave(1).z;
    nfreq = numel(freq);
    period = 1./freq;

    % Impedances    
    Zxx = cat(1,ffsave.Zxx);                Zxx = smooth(log10(period),Zxx(ind),'loess');
    Zxx_Err = cat(1,ffsave.Zxx_Err);        Zxx_Err = Zxx_Err(ind);
    Zxy = cat(1,ffsave.Zxy);                Zxy = smooth(log10(period),Zxy(ind),'loess');
    Zxy_Err = cat(1,ffsave.Zxy_Err);        Zxy_Err = Zxy_Err(ind);
    Zyx = cat(1,ffsave.Zyx);                Zyx = smooth(log10(period),Zyx(ind),'loess');
    Zyx_Err = cat(1,ffsave.Zyx_Err);        Zyx_Err = Zyx_Err(ind);
    Zyy = cat(1,ffsave.Zyy);                Zyy = smooth(log10(period),Zyy(ind),'loess');
    Zyy_Err = cat(1,ffsave.Zyy_Err);        Zyy_Err = Zyy_Err(ind);

%     Zxx = cat(1,ffsave.Zxx);                Zxx = Zxx(ind);
%     Zxx_Err = cat(1,ffsave.Zxx_Err);        Zxx_Err = Zxx_Err(ind);
%     Zxy = cat(1,ffsave.Zxy);                Zxy = Zxy(ind);
%     Zxy_Err = cat(1,ffsave.Zxy_Err);        Zxy_Err = Zxy_Err(ind);
%     Zyx = cat(1,ffsave.Zyx);                Zyx = Zyx(ind);
%     Zyx_Err = cat(1,ffsave.Zyx_Err);        Zyx_Err = Zyx_Err(ind);
%     Zyy = cat(1,ffsave.Zyy);                Zyy = Zyy(ind);
%     Zyy_Err = cat(1,ffsave.Zyy_Err);        Zyy_Err = Zyy_Err(ind);

    % Resistivity, Phase and Errors
    [rho,phi] = rhoandphi(Zxx,Zxy,Zyx,Zyy,Zxx_Err,Zxy_Err,Zyx_Err,Zyy_Err,freq);

    % Tipper
    txz = cat(1,ffsave.txz);                txz = smooth(log10(period),txz(ind),'loess');
    txz_Err = cat(1,ffsave.txz_Err);        txz_Err = txz_Err(ind);
    tyz = cat(1,ffsave.tyz);                tyz = smooth(log10(period),tyz(ind),'loess');    
    tyz_Err = cat(1,ffsave.tyz_Err);        tyz_Err = tyz_Err(ind);

%     txz = cat(1,ffsave.txz);                txz = txz(ind);
%     txz_Err = cat(1,ffsave.txz_Err);        txz_Err = txz_Err(ind);
%     tyz = cat(1,ffsave.tyz);                tyz = tyz(ind);    
%     tyz_Err = cat(1,ffsave.tyz_Err);        tyz_Err = tyz_Err(ind);

    % HMTF
    try
        txx = cat(1,ffsave.txx);                txx = txx(ind);
        txy = cat(1,ffsave.txy);                txy = txy(ind);
        tyx = cat(1,ffsave.tyx);                tyx = tyx(ind);
        tyy = cat(1,ffsave.tyy);                tyy = tyy(ind);
        t_id = cat(1,ffsave.t_id);              t_id = t_id(ind,:);
    catch
    end

    % Saving site variable FFsave
    FFsave = create_mt;
    % Metainfo
    FFsave.source = {'data'};
    FFsave.lonlat = lonlat;
    FFsave.UTM = [x,y];
    FFsave.z = z;
    FFsave.nfreq = nfreq;
    FFsave.freq = freq;
    FFsave.per = period;
    % Impedances
    FFsave.Zxx = Zxx;                       FFsave.Zxx_Err = Zxx_Err;
    FFsave.Zxy = Zxy;                       FFsave.Zxy_Err = Zxy_Err;
    FFsave.Zyx = Zyx;                       FFsave.Zyx_Err = Zyx_Err;
    FFsave.Zyy = Zyy;                       FFsave.Zyy_Err = Zyy_Err;
    % Resistivities
    FFsave.rhoxx = rho.rhoxx;               FFsave.rhoxy = rho.rhoxy;    
    FFsave.rhoyx = rho.rhoyx;               FFsave.rhoyy = rho.rhoyy;
    FFsave.rhoxx_Err = rho.rhoxx_Err;       FFsave.rhoxy_Err = rho.rhoxy_Err;    
    FFsave.rhoyx_Err = rho.rhoyx_Err;       FFsave.rhoyy_Err = rho.rhoyy_Err;
    % Phases
    FFsave.phixx = phi.phixx;               FFsave.phixy = phi.phixy;    
    FFsave.phiyx = phi.phiyx;               FFsave.phiyy = phi.phiyy;
    FFsave.phixx_Err = phi.phixx_Err;       FFsave.phixy_Err = phi.phixy_Err;    
    FFsave.phiyx_Err = phi.phiyx_Err;       FFsave.phiyy_Err = phi.phiyy_Err;
    % Tipper
    FFsave.txz = txz;                       FFsave.tyz = tyz;
    FFsave.txz_Err = txz_Err;               FFsave.tyz_Err = tyz_Err;

    % Calculating Tensors
    FFsave = CART(FFsave);
    % Correcting Phase Tensors
    FFsave = calc_PT(FFsave);
    
    % Saving file
    if max(freq) >= 1
        f1 = round(max(freq));
        unit1 = 'Hz';
    else
        f1 = round(1./max(freq));
        unit1 = 's';
    end
    if min(freq) >= 1
        f2 = round(min(freq));
        unit2 = 'Hz';
    else
        f2 = round(1./min(freq));
        unit2 = 's';
    end
    range = [num2str(f1),unit1,'-',num2str(f2),unit2];

    if rr
        RR = 'RR';
    else
        RR = [];
    end

    % Filename  
    filename = [name{:},RR,'_',range,'_',model{:},'_',reg{:},'_',type{:},'.mat'];
    FFsave.site = filename;
    mt = FFsave;
    save(filename,'mt');
    
end

%% Function to create RAW files
function mergeRAW(filelist,model,reg,rr)
     for i = 1:numel(filelist)
        load(filelist{i})
        job_merge{i} = job;
        act_chan_merge{i} = act_chan; 
        RAW_merge{i} = RAW;
        EV_merge{i} = EV;        
     end

    for i = 1:numel(EV_merge)
        sr = EV_merge{i}.sr;
    end
    [~,sr_id] = sort(sr,'descend');
    
    % Arranging files
    RAW_merge = RAW_merge(sr_id);
    EV_merge = EV_merge(sr_id);
    
    % Ammount of best ranked EV
    for i = 1:numel(EV_merge)
        perc = EV_merge{1}.perc;
    end
    
    if numel(unique(perc)) ~= 1
        perc = min(perc);
    end
    
    % Cleaning memory
    clear job act_chan RAW 
    
    % Concatenating RAW structure
    RAW = create_ffsave;
    for j = 1:perc
        
        % Allocating variables
        site = [];      lonlat = [];    UTM = [];       z = [];
        nfreq = [];     freq = [];      per = [];
        Zxx = [];       Zxy = [];       Zyx = [];       Zyy = [];
        Zxx_Err = [];   Zxy_Err = [];   Zyx_Err = [];   Zyy_Err = [];
        txz = [];       txz_Err = [];
        tyz = [];       tyz_Err = [];   
        txx = [];       txy = [];       tyx = [];       tyy = [];
        t_id = [];
        phi11 = [];     phi12 = [];     phi21 = [];     phi22 = [];        
        alpha = [];     beta = [];      theta = [];     lambda = [];
        alpha_Err = []; beta_Err = [];  theta_Err = []; lambda_Err = [];
        phimax = [];    phimin = [];        
        phimax_Err = [];phimin_Err = [];        
        
        for i = 1:numel(RAW_merge)
            
            % Name and coordinates
            site = [site;RAW_merge{i}(j).site];
            lonlat = [lonlat;{RAW_merge{i}(j).lonlat}];
            UTM = [UTM;{RAW_merge{i}(j).UTM}];
            z = [z;{RAW_merge{i}(j).z}];
            
            % Frequencies
            nfreq = [nfreq;RAW_merge{i}(j).nfreq];
            freq = [freq;RAW_merge{i}(j).freq];
            per = [per;RAW_merge{i}(j).per];
            
            % Impedances
            Zxx = [Zxx;RAW_merge{i}(j).Zxx];
            Zxy = [Zxy;RAW_merge{i}(j).Zxy];
            Zyx = [Zyx;RAW_merge{i}(j).Zyx];
            Zyy = [Zyy;RAW_merge{i}(j).Zyy];
            Zxx_Err = [Zxx_Err;RAW_merge{i}(j).Zxx_Err];
            Zxy_Err = [Zxy_Err;RAW_merge{i}(j).Zxy_Err];
            Zyx_Err = [Zyx_Err;RAW_merge{i}(j).Zyx_Err];
            Zyy_Err = [Zyy_Err;RAW_merge{i}(j).Zyy_Err];
            
            % Tipper
            txz = [txz;RAW_merge{i}(j).txz];
            tyz = [tyz;RAW_merge{i}(j).tyz];
            txz_Err = [txz_Err;RAW_merge{i}(j).txz_Err];
            tyz_Err = [tyz_Err;RAW_merge{i}(j).tyz_Err];
            
            % Phase Tensor
            phi11 = [phi11;RAW_merge{i}(j).phi11];
            phi12 = [phi12;RAW_merge{i}(j).phi12];
            phi21 = [phi21;RAW_merge{i}(j).phi21];
            phi22 = [phi22;RAW_merge{i}(j).phi22];
            alpha = [alpha;RAW_merge{i}(j).alpha];
            beta = [beta;RAW_merge{i}(j).beta];
            theta = [theta;RAW_merge{i}(j).theta];
            lambda = [lambda;RAW_merge{i}(j).lambda];
            phimax = [phimax;RAW_merge{i}(j).phimax];
            phimin = [phimin;RAW_merge{i}(j).phimin];
            alpha_Err = [alpha_Err;RAW_merge{i}(j).alpha_Err];
            beta_Err = [beta_Err;RAW_merge{i}(j).beta_Err];
            theta_Err = [theta_Err;RAW_merge{i}(j).theta_Err];
            lambda_Err = [lambda_Err;RAW_merge{i}(j).lambda_Err];
            phimax_Err = [phimax_Err;RAW_merge{i}(j).phimax_Err];
            phimin_Err = [phimin_Err;RAW_merge{i}(j).phimin_Err];
            
            % HMTF
            txx = [txx;RAW_merge{i}(j).txx];       
            txy = [txy;RAW_merge{i}(j).txy]; 
            tyx = [tyx;RAW_merge{i}(j).tyx];         
            tyy = [tyy;RAW_merge{i}(j).tyy]; 
            t_id = [t_id;RAW_merge{i}(j).t_id];
            
        end
        
        % Name and coordinates
        RAW(j).site = unique(site);
        RAW(j).nfreq = unique(nfreq);
        RAW(j).lonlat = lonlat{1};
        RAW(j).UTM = UTM{1};
        RAW(j).z = z{1};
        
        % Frequencies
        RAW(j).nfreq = sum(nfreq);
        RAW(j).freq = freq;
        RAW(j).per = per;
        
        % Saving transfer functions in RAW structure: Z
        RAW(j).Zxx = Zxx;                           RAW(j).Zxx_Err = Zxx_Err;                       
        RAW(j).Zxy = Zxy;                           RAW(j).Zxy_Err = Zxy_Err;
        RAW(j).Zyx = Zyx;                           RAW(j).Zyx_Err = Zyx_Err;        
        RAW(j).Zyy = Zyy;                           RAW(j).Zyy_Err = Zyy_Err;   
        
        % Saving transfer functions in RAW structure: T
        RAW(j).txz = txz;                            RAW(j).txz_Err = txz_Err;
        RAW(j).tyz = tyz;                            RAW(j).tyz_Err = tyz_Err;
        
        % Saving transfer functions in RAW structure: PT
        RAW(j).phi11 = phi11;                        RAW(j).phi12 = phi12;
        RAW(j).phi21 = phi21;                        RAW(j).phi22 = phi22;
        RAW(j).beta = beta;                          RAW(j).beta_Err = beta_Err;
        RAW(j).alpha = alpha;                        RAW(j).alpha_Err = alpha_Err;
        RAW(j).theta = theta;                        RAW(j).theta_Err = theta_Err;
        RAW(j).lambda = lambda;                      RAW(j).lambda_Err = lambda_Err;
        RAW(j).phimax = phimax;                      RAW(j).phimax_Err = phimax_Err;
        RAW(j).phimin = phimin;                      RAW(j).phimin_Err = phimin_Err;   
        
    end
    
    % Resistivity and phase calculations, scaling and corrections
    RAW = FFproc_RAW_rhoandphi(RAW,job_merge{1});   
    
    clear EV
    % Concatenating EV structure
    for j = 2:perc
        best = [];      
        for i = 1:numel(RAW_merge)
            % Best estimations
            best = [best;EV_merge{i}.best(j).best];
        end
        EV.best(j).best = best;
    end
    
    % Allocating variables
    perc = [];      EVal = [];      Ev = []; 
    timeseg = [];   tfreq = [];     ntseg = []; 
    tstart = [];    tend = [];  
    
    for i = 1:numel(RAW_merge)
        
        perc = [perc,EV_merge{i}.perc];
        EVal = [EVal,EV_merge{i}.EVal]; 
        Ev = [Ev,EV_merge{i}.EV]; 
        timeseg = [timeseg,EV_merge{i}.timeseg]; 
        tfreq = [tfreq;EV_merge{i}.tfreq]; 
        ntseg = [ntseg;EV_merge{i}.tfreq];
        tstart = [tstart;EV_merge{i}.tstart];
        tend = [tend;EV_merge{i}.tend];
        
    end
    
    EV.perc = min(perc);
    EV.EVal = EVal;
    EV.EV = EV;
    EV.timeseg = timeseg;
    EV.tfreq = tfreq;
    EV.ntseg = ntseg;
    EV.tstart = tstart;
    EV.tend = tend;      
    
    % Saving magnetic reference status (to identify local or RR)
    for i = 1:numel(job_merge)
        for j = 1:numel(job_merge{i}.ffts)
            magref{i,j} = [job_merge{i}.ffts(j).rsbx,job_merge{i}.ffts(j).rsby];
            device{i,j} = job_merge{i}.ffts(j).device;
        end
    end  
    
    site = unique(site);    
    for i = 1:size(magref,2)
        mag(i,:) = all(cell2mat(magref(:,i)));     
        job.ffts.rsbx = mag(i,1);
        job.ffts.rsby = mag(i,2);
        job.name = site;
        job.ffts.GPS = [lonlat(i,:),z];
        job.ffts.device = char(unique(device));
        job.eigperc = min(perc);
        job.ftarg = 6;
        job.ffts.t = [];
        job.ffts.data = [];
    end
    
    % Saving active channels array (to identify which channels were used
    % to perform calculations
    for i = 1:numel(job_merge)
        for j = 1:numel(job_merge{i}.ffts)
            act{i,j} = act_chan_merge{1,i};
        end
    end
    for i = 1:size(act,2)
        act_chan{i} = all(reshape(cell2mat([act{:,i}]),5,[])');
    end    
    
    % Retrieving target frequency range label: seconds or Hz
    if max(freq) >= 1
        f1 = round(max(freq));
        unit1 = 'Hz';
    else
        f1 = round(1./max(freq));
        unit1 = 's';
    end
    if min(freq) >= 1
        f2 = round(min(freq));
        unit2 = 'Hz';
    else
        f2 = round(1./min(freq));
        unit2 = 's';
    end

    switch job_merge{1}.model
        case 'Standard'
            model = 'STD_';
        case 'Advanced'
            model = 'ANM_';
    end

    switch job_merge{1}.reg
        case 'Robust'
            reg = 'ROB_';
        case 'LSQ'
            reg = 'LSQ_';
    end

    if rr
        RR = '_RR';
    else
        RR = [];
    end
    
    % Saving file
    outname = ['AUTORAW_',model,reg,num2str(f1),unit1,'-',num2str(f2),unit2,RR,'_RAW',num2str(min(perc)),'%'];       
    Results = RAW;
    target_freq = tfreq;
    save(outname,'Results','EV','job','act_chan','target_freq');
end