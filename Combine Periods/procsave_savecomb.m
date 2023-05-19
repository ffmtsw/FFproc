function [filename,path] = procsave_savecomb(app,ffsave,ev)
    
    % Indices to be included
    include = find(cat(1,app.table.Data{:,3}));

    siteid = [];
    siteind = [];
    for  i = 1:numel(ffsave)    
        siteid = [siteid,ones(1,ffsave(i).nfreq).*i];
        siteind = [siteind,1:ffsave(i).nfreq];
    end
    idcheck = [siteid',siteind'];
    ID = str2num([app.table.Data{:,1}]'); 
    ID = ID(include);

    perstr = app.table.Data(:,2);
    % Selecting the right label (scale Hz or s)
    per = zeros(numel(perstr),1);
    for i = 1:numel(perstr)
        if strcmp(perstr{i}(end-3:end),'(Hz)')
            per(i) = 1/str2double(perstr{i}(1:end-4));
        elseif strcmp(perstr{i}(end-2:end),'(s)')
            per(i) = str2double(perstr{i}(1:end-3));
        end
    end
    
    if numel(sort(round(per(include),4))) ~= numel(unique(round(per,4))) |...
             sort(round(per(include),4)) ~= unique(round(per,4))
         
        answer = uiconfirm(app.UIFigure,{'Site contains repeated periods';...
                                         'Coincident data will be replaced for those estimated at longer periods'},...
                                         'FFproc - Combine Periods',...
                                         'Options',{'Continue','Return'},'DefaultOption',1,'CancelOption',2,...
                                         'Icon','Warning');
        switch answer
            case 'Continue'
                [period,ind,~] = unique(round(per,4),'stable');
                % Metainfo
                lonlat = mean(cat(1,ffsave.lonlat));
                [x,y,~] = deg2utm(lonlat(1),lonlat(2));
                z = mean(cat(1,ffsave.z));
                nfreq = numel(period);
                freq = 1./period;
                % Impedances
                Zxx = cat(1,ffsave.Zxx);                Zxx = Zxx(ind);
                Zxx_Err = cat(1,ffsave.Zxx_Err);        Zxx_Err = Zxx_Err(ind);
                Zxy = cat(1,ffsave.Zxy);                Zxy = Zxy(ind);
                Zxy_Err = cat(1,ffsave.Zxy_Err);        Zxy_Err = Zxy_Err(ind);
                Zyx = cat(1,ffsave.Zyx);                Zyx = Zyx(ind);
                Zyx_Err = cat(1,ffsave.Zyx_Err);        Zyx_Err = Zyx_Err(ind);
                Zyy = cat(1,ffsave.Zyy);                Zyy = Zyy(ind);
                Zyy_Err = cat(1,ffsave.Zyy_Err);        Zyy_Err = Zyy_Err(ind);
                % Resistivities
                rhoxx = cat(1,ffsave.rhoxx);            rhoxx = rhoxx(ind);
                rhoxy = cat(1,ffsave.rhoxy);            rhoxy = rhoxy(ind);
                rhoyx = cat(1,ffsave.rhoyx);            rhoyx = rhoyx(ind);
                rhoyy = cat(1,ffsave.rhoyy);            rhoyy = rhoyy(ind);
                rhoxx_Err = cat(1,ffsave.rhoxx_Err);    rhoxx_Err = rhoxx_Err(ind);
                rhoxy_Err = cat(1,ffsave.rhoxy_Err);	rhoxy_Err = rhoxy_Err(ind);
                rhoyx_Err = cat(1,ffsave.rhoyx_Err);	rhoyx_Err = rhoyx_Err(ind);
                rhoyy_Err = cat(1,ffsave.rhoyy_Err);	rhoyy_Err = rhoyy_Err(ind);
                % Phases
                phixx = cat(1,ffsave.phixx);            phixx = phixx(ind);
                phixy = cat(1,ffsave.phixy);            phixy = phixy(ind);
                phiyx = cat(1,ffsave.phiyx);            phiyx = phiyx(ind);
                phiyy = cat(1,ffsave.phiyy);            phiyy = phiyy(ind);
                phixx_Err = cat(1,ffsave.phixx_Err);	phixx_Err = phixx_Err(ind);
                phixy_Err = cat(1,ffsave.phixy_Err);	phixy_Err = phixy_Err(ind);
                phiyx_Err = cat(1,ffsave.phiyx_Err);	phiyx_Err = phiyx_Err(ind);
                phiyy_Err = cat(1,ffsave.phiyy_Err);	phiyy_Err = phiyy_Err(ind);
                % Tipper
                txz = cat(1,ffsave.txz);                txz = txz(ind);
                tyz = cat(1,ffsave.tyz);                tyz = tyz(ind);
                txz_Err = cat(1,ffsave.txz_Err);        txz_Err = txz_Err(ind);
                tyz_Err = cat(1,ffsave.tyz_Err);        tyz_Err = tyz_Err(ind);
                % PT
                phi11 = cat(1,ffsave.phi11);            phi11 = phi11(ind);
                phi12 = cat(1,ffsave.phi12);            phi12 = phi12(ind);
                phi21 = cat(1,ffsave.phi21);            phi21 = phi21(ind);
                phi22 = cat(1,ffsave.phi22);            phi22 = phi22(ind);
                beta = cat(1,ffsave.beta);              beta = beta(ind);
                beta_Err = cat(1,ffsave.beta_Err);      beta_Err = beta_Err(ind);
                alpha = cat(1,ffsave.alpha);            alpha = alpha(ind);
                alpha_Err = cat(1,ffsave.alpha_Err);	alpha_Err = alpha_Err(ind);
                theta = cat(1,ffsave.theta);            theta = theta(ind);
                theta_Err = cat(1,ffsave.theta_Err);	theta_Err = theta_Err(ind);
                lambda = cat(1,ffsave.lambda);          lambda = lambda(ind);
                lambda_Err = cat(1,ffsave.lambda_Err);	lambda_Err = lambda_Err(ind);
                phimax = cat(1,ffsave.phimax);          phimax = phimax(ind);
                phimax_Err = cat(1,ffsave.phimax_Err);	phimax_Err = phimax_Err(ind);
                phimin = cat(1,ffsave.phimin);          phimin = phimin(ind);
                phimin_Err = cat(1,ffsave.phimin_Err);	phimin_Err = phimin_Err(ind);
                % HMTF
                txx = cat(1,ffsave.txx); 
                txy = cat(1,ffsave.txy);
                tyx = cat(1,ffsave.tyx);
                tyy = cat(1,ffsave.tyy); 
                t_id = cat(1,ffsave.t_id);
                if ~isempty(txx)
                    txx = txx(ind);
                    txy = txy(ind);
                    tyx = tyx(ind);
                    tyy = tyy(ind);
                    t_id = t_id(ind,:);
                end
                idcheck = idcheck(ind,:);
            case 'Return'
                return
        end
    else
        [period,ind] = unique(per(include),'stable');
        % Metainfo
        coors = cat(1,ffsave.lonlat);
        lonlat = [mean(coors(:,1)),mean(coors(:,2))];
        [x,y,~] = deg2utm(lonlat(1),lonlat(2));
        z = mean(cat(1,ffsave.z));
        nfreq = numel(period);
        freq = 1./period;
        % Impedances
        Zxx = cat(1,ffsave.Zxx);                Zxx = Zxx(include);                 Zxx = Zxx(ind);
        Zxx_Err = cat(1,ffsave.Zxx_Err);        Zxx_Err = Zxx_Err(include);         Zxx_Err = Zxx_Err(ind);
        Zxy = cat(1,ffsave.Zxy);                Zxy = Zxy(include);                 Zxy = Zxy(ind);
        Zxy_Err = cat(1,ffsave.Zxy_Err);        Zxy_Err = Zxy_Err(include);         Zxy_Err = Zxy_Err(ind);
        Zyx = cat(1,ffsave.Zyx);                Zyx = Zyx(include);                 Zyx = Zyx(ind);
        Zyx_Err = cat(1,ffsave.Zyx_Err);        Zyx_Err = Zyx_Err(include);         Zyx_Err = Zyx_Err(ind);
        Zyy = cat(1,ffsave.Zyy);                Zyy = Zyy(include);                 Zyy = Zyy(ind);
        Zyy_Err = cat(1,ffsave.Zyy_Err);        Zyy_Err = Zyy_Err(include);         Zyy_Err = Zyy_Err(ind);
        % Resistivities
        rhoxx = cat(1,ffsave.rhoxx);            rhoxx = rhoxx(include);             rhoxx = rhoxx(ind);
        rhoxy = cat(1,ffsave.rhoxy);            rhoxy = rhoxy(include);             rhoxy = rhoxy(ind);
        rhoyx = cat(1,ffsave.rhoyx);            rhoyx = rhoyx(include);             rhoyx = rhoyx(ind);
        rhoyy = cat(1,ffsave.rhoyy);            rhoyy = rhoyy(include);             rhoyy = rhoyy(ind);
        rhoxx_Err = cat(1,ffsave.rhoxx_Err);    rhoxx_Err = rhoxx_Err(include);     rhoxx_Err = rhoxx_Err(ind);
        rhoxy_Err = cat(1,ffsave.rhoxy_Err);    rhoxy_Err = rhoxy_Err(include);     rhoxy_Err = rhoxy_Err(ind);
        rhoyx_Err = cat(1,ffsave.rhoyx_Err);    rhoyx_Err = rhoyx_Err(include);     rhoyx_Err = rhoyx_Err(ind);
        rhoyy_Err = cat(1,ffsave.rhoyy_Err);    rhoyy_Err = rhoyy_Err(include);     rhoyy_Err = rhoyy_Err(ind);
        % Phases
        phixx = cat(1,ffsave.phixx);            phixx = phixx(include);             phixx = phixx(ind);
        phixy = cat(1,ffsave.phixy);            phixy = phixy(include);             phixy = phixy(ind);
        phiyx = cat(1,ffsave.phiyx);            phiyx = phiyx(include);             phiyx = phiyx(ind);
        phiyy = cat(1,ffsave.phiyy);            phiyy = phiyy(include);             phiyy = phiyy(ind);
        phixx_Err = cat(1,ffsave.phixx_Err);    phixx_Err = phixx_Err(include);     phixx_Err = phixx_Err(ind);
        phixy_Err = cat(1,ffsave.phixy_Err);    phixy_Err = phixy_Err(include);     phixy_Err = phixy_Err(ind);
        phiyx_Err = cat(1,ffsave.phiyx_Err);    phiyx_Err = phiyx_Err(include);     phiyx_Err = phiyx_Err(ind);
        phiyy_Err = cat(1,ffsave.phiyy_Err);    phiyy_Err = phiyy_Err(include);     phiyy_Err = phiyy_Err(ind);
        % Tipper
        txz = cat(1,ffsave.txz);                txz = txz(include);                 txz = txz(ind);
        tyz = cat(1,ffsave.tyz);                tyz = tyz(include);                 tyz = tyz(ind);
        txz_Err = cat(1,ffsave.txz_Err);        txz_Err = txz_Err(include);         txz_Err = txz_Err(ind);
        tyz_Err = cat(1,ffsave.tyz_Err);        tyz_Err = tyz_Err(include);         tyz_Err = tyz_Err(ind);
        % PT
        phi11 = cat(1,ffsave.phi11);            phi11 = phi11(include);             phi11 = phi11(ind);
        phi12 = cat(1,ffsave.phi12);            phi12 = phi12(include);             phi12 = phi12(ind);
        phi21 = cat(1,ffsave.phi21);            phi21 = phi21(include);             phi21 = phi21(ind);
        phi22 = cat(1,ffsave.phi22);            phi22 = phi22(include);             phi22 = phi22(ind);
        beta = cat(1,ffsave.beta);              beta = beta(include);               beta = beta(ind);
        beta_Err = cat(1,ffsave.beta_Err);      beta_Err = beta_Err(include);       beta_Err = beta_Err(ind);
        alpha = cat(1,ffsave.alpha);            alpha = alpha(include);             alpha = alpha(ind);
        alpha_Err = cat(1,ffsave.alpha_Err);    alpha_Err = alpha_Err(include);     alpha_Err = alpha_Err(ind);
        theta = cat(1,ffsave.theta);            theta = theta(include);             theta = theta(ind);
        theta_Err = cat(1,ffsave.theta_Err);    theta_Err = theta_Err(include);     theta_Err = theta_Err(ind);
        lambda = cat(1,ffsave.lambda);          lambda = lambda(include);           lambda = lambda(ind);
        lambda_Err = cat(1,ffsave.lambda_Err);  lambda_Err = lambda_Err(include);   lambda_Err = lambda_Err(ind);
        phimax = cat(1,ffsave.phimax);          phimax = phimax(include);           phimax = phimax(ind);
        phimax_Err = cat(1,ffsave.phimax_Err);  phimax_Err = phimax_Err(include);   phimax_Err = phimax_Err(ind);
        phimin = cat(1,ffsave.phimin);          phimin = phimin(include);           phimin = phimin(ind);
        phimin_Err = cat(1,ffsave.phimin_Err);  phimin_Err = phimin_Err(include);   phimin_Err = phimin_Err(ind);
        try
            % HMTF
            txx = cat(1,ffsave.txx);                txx = txx(include);                 txx = txx(ind);                        
            txy = cat(1,ffsave.txy);                txy = txy(include);                 txy = txy(ind);
            tyx = cat(1,ffsave.tyx);                tyx = tyx(include);                 tyx = tyx(ind);
            tyy = cat(1,ffsave.tyy);                tyy = tyy(include);                 tyy = tyy(ind);
            t_id = cat(1,ffsave.t_id);              t_id = t_id(include,:);             t_id = t_id(ind,:);
        catch
        end
        idcheck = idcheck(ind,:);
        
    end
    
    %% Saving site variable FFsave
    FFsave = create_ffsave;
    % Metainfo
    FFsave.lonlat = lonlat;
    FFsave.UTM = [x,y];
    FFsave.z = z;
    FFsave.nfreq = nfreq;
    FFsave.freq = freq;
    FFsave.per = period;
    % Impedances
    FFsave.Zxx = Zxx;               FFsave.Zxx_Err = Zxx_Err;
    FFsave.Zxy = Zxy;               FFsave.Zxy_Err = Zxy_Err;
    FFsave.Zyx = Zyx;               FFsave.Zyx_Err = Zyx_Err;
    FFsave.Zyy = Zyy;               FFsave.Zyy_Err = Zyy_Err;
    % Resistivities
    FFsave.rhoxx = rhoxx;           FFsave.rhoxy = rhoxy;    
    FFsave.rhoyx = rhoyx;           FFsave.rhoyy = rhoyy;
    FFsave.rhoxx_Err = rhoxx_Err;   FFsave.rhoxy_Err = rhoxy_Err;    
    FFsave.rhoyx_Err = rhoyx_Err;   FFsave.rhoyy_Err = rhoyy_Err;
    % Phases
    FFsave.phixx = phixx;           FFsave.phixy = phixy;    
    FFsave.phiyx = phiyx;           FFsave.phiyy = phiyy;
    FFsave.phixx_Err = phixx_Err;   FFsave.phixy_Err = phixy_Err;    
    FFsave.phiyx_Err = phiyx_Err;   FFsave.phiyy_Err = phiyy_Err;
    % Tipper
    FFsave.txz = txz;               FFsave.tyz = tyz;
    FFsave.txz_Err = txz_Err;       FFsave.tyz_Err = tyz_Err;
    % PT
    FFsave.phi11 = phi11;           FFsave.phi12 = phi12;
    FFsave.phi21 = phi21;           FFsave.phi22 = phi22;
    FFsave.phimin = phimin;         FFsave.phimin_Err = phimin_Err;
    FFsave.phimax = phimax;         FFsave.phimax_Err = phimax_Err;     
    FFsave.alpha = alpha;           FFsave.alpha_Err = alpha_Err;   
    FFsave.beta = beta;             FFsave.beta_Err = beta_Err;
    FFsave.theta = theta;           FFsave.theta_Err = theta_Err;
    FFsave.lambda = lambda;         FFsave.lambda_Err = lambda_Err;
    try
        % HMTF
        FFsave.txx = txx;
        FFsave.txy = txy;
        FFsave.tyx = tyx;
        FFsave.tyy = tyy;
        FFsave.t_id = t_id;
    catch
    end
    
    %% Saving Eigenvalue variable EV
    try
        n = 1;
        m = 0;
        for i = 1:numel(idcheck(:,1))
            if n == idcheck(i,1)
                m = m+1;
            else
                m = 1;
                n = n+1;
            end
            EV(n).best(m,:) = ev(idcheck(i,1)).best(idcheck(i,2),:);
            EV(n).EVal(m) = ev(idcheck(i,1)).EVal(idcheck(i,2));
            EV(n).EV(m) = ev(idcheck(i,1)).EV(idcheck(i,2));
            EV(n).timeseg(m) = ev(idcheck(i,1)).timeseg(idcheck(i,2));
            EV(n).tfreq(m) = ev(idcheck(i,1)).tfreq(idcheck(i,2));
            EV(n).ntseg(m) = ev(idcheck(i,1)).ntseg(idcheck(i,2));
            EV(n).tseg(m) = ev(idcheck(i,1)).tseg(idcheck(i,2));
        end
        for i = 1:numel(ev)
            EV(i).perc = ev(i).perc;
            EV(i).sr = ev(i).sr;
            EV(i).tstart = ev(i).tstart;
            EV(i).tend = ev(i).tend;
        end    
        EV = orderfields(EV,ev);
    catch
        EV = [];
    end
       
    %% Saving file
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
    filename = unique(app.jobtable.Data(:,2));
    filename = [[filename{:}],'_',range];
    
    [filename,path] = uiputfile(['COMB_',filename,'.mat'],'Save site variable with combined periods in *.mat file');
    if path == 0
        return
    else
        cd(path)    
        FFsave.site = {filename};        
        save(filename,'FFsave','EV');
    end    
    
end