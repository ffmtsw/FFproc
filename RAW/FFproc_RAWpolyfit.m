function FFproc_RAWpolyfit(job,RAW,EV,act_chan,eigperc,polyord)
    
    nts = numel(job.ffts);
    target_freq = EV.tfreq;
    
    % Magnetic reference channels
    for i = 1:numel(job.ffts)
        magref(i,1) = job.ffts(i).rsbx;
        magref(i,2) = job.ffts(i).rsby;
    end
    
    % Create FFsave structure
    FFsave = ffsave_create;

    for a = 1:nts % Number of time series per job
        segments = repmat(act_chan,eigperc-1,1).';
%         if act_chan{a}(1) && act_chan{a}(2)
        if act_chan{a}(3) && act_chan{a}(4)
            %% Extracting info from Results structure variable
            for l = 1:numel(target_freq)  % Number of time segments per job
                for w = 2:eigperc  % Eigenvalue percentage
                    zxx(l,w-1) = RAW(w).Zxx(l,a);
                    zxy(l,w-1) = RAW(w).Zxy(l,a);
                    zyx(l,w-1) = RAW(w).Zyx(l,a);
                    zyy(l,w-1) = RAW(w).Zyy(l,a);
                    zxx_Err(l,w-1) = RAW(w).Zxx_Err(l,a);
                    zxy_Err(l,w-1) = RAW(w).Zxy_Err(l,a);
                    zyx_Err(l,w-1) = RAW(w).Zyx_Err(l,a);
                    zyy_Err(l,w-1) = RAW(w).Zyy_Err(l,a);
                end
            end
            %% Polynom fitting
            % Zxx
            [~,comp,err,~,~,~] = fitpoly2(1./target_freq,real(zxx),eigperc,real(zxx_Err),0,1,polyord);
            Zxx_Re = comp;  % Real
            FFsave(a).Zxx_Re = Zxx_Re';
            [~,comp,~,~,~,~] = fitpoly2(1./target_freq,imag(zxx),eigperc,imag(zxx_Err),0,1,polyord);
            Zxx_Im = comp;  % Imaginary
            FFsave(a).Zxx_Im = Zxx_Im';
            FFsave(a).Zxx_Err = err';         % Error
            FFsave(a).Zxx = complex(FFsave(a).Zxx_Re,FFsave(a).Zxx_Im);  % Complex
            
            % Zxy
            m = median(sign(real(zxy(:))));
            [~,comp,err,~,~,~] = fitpoly2(1./target_freq,log10(abs(real(zxy))),eigperc,real(zxy_Err),0,1,polyord);
            Zxy_Re = m.*10.^comp;  % Real
            FFsave(a).Zxy_Re = Zxy_Re';
            m = median(sign(imag(zxy(:))));
            [~,comp,~,~,~,~] = fitpoly2(1./target_freq,log10(abs(imag(zxy))),eigperc,imag(zxy_Err),0,1,polyord);
            Zxy_Im = m.*10.^comp;  % Imaginary
            FFsave(a).Zxy_Im = Zxy_Im';
            FFsave(a).Zxy_Err = err';         % Error
            FFsave(a).Zxy = complex(FFsave(a).Zxy_Re,FFsave(a).Zxy_Im);  % Complex
            
            % Zyx
            m = median(sign(real(zyx(:))));
            [~,comp,err,~,~,~] = fitpoly2(1./target_freq,log10(abs(real(zyx))),eigperc,real(zyx_Err),0,1,polyord);
            Zyx_Re = m.*10.^comp;  % Real
            FFsave(a).Zyx_Re = Zyx_Re';
            m = median(sign(imag(zyx(:))));
            [~,comp,~,~,~,~] = fitpoly2(1./target_freq,log10(abs(imag(zyx))),eigperc,imag(zyx_Err),0,1,polyord);
            Zyx_Im = m.*10.^comp;  % Imaginary
            FFsave(a).Zyx_Im = Zyx_Im';
            FFsave(a).Zyx_Err = err';         % Error
            FFsave(a).Zyx = complex(FFsave(a).Zyx_Re,FFsave(a).Zyx_Im);  % Complex
            
            % Zyy
            [~,comp,err,~,~,~] = fitpoly2(1./target_freq,real(zyy),eigperc,real(zyy_Err),0,1,polyord);
            Zyy_Re = comp;  % Real
            FFsave(a).Zyy_Re = Zyy_Re';
            [~,comp,~,~,~,~] = fitpoly2(1./target_freq,imag(zyy),eigperc,imag(zyy_Err),0,1,polyord);
            Zyy_Im = comp;  % Imaginary
            FFsave(a).Zyy_Im = Zyy_Im';
            FFsave(a).Zyy_Err = err';         % Error
            FFsave(a).Zyy = complex(FFsave(a).Zyy_Re,FFsave(a).Zyy_Im);  % Complex
            
            %% Phase Tensor
            for l = 1:numel(target_freq)
                Z = [FFsave(a).Zxx(l);FFsave(a).Zyx(l);FFsave(a).Zxy(l);FFsave(a).Zyy(l)];
                Z_err = [FFsave(a).Zxx_Err(l);FFsave(a).Zyx_Err(l);...
                         FFsave(a).Zxy_Err(l);FFsave(a).Zyy_Err(l)];
                PT = FFproc_PT(Z);
                PT_err = FFproc_PT_err(Z,Z_err);
                FFsave(a).phi11(l,1) = PT.phi11;        FFsave(a).phi12(l,1) = PT.phi12;
                FFsave(a).phi21(l,1) = PT.phi21;        FFsave(a).phi22(l,1) = PT.phi22;
                FFsave(a).phimin(l,1) = PT.phimin;      FFsave(a).phimax(l,1) = PT.phimax;
                FFsave(a).alpha(l,1) = PT.alpha;        FFsave(a).beta(l,1) = PT.beta;
                FFsave(a).theta(l,1) = PT.theta;        FFsave(a).lambda(l,1) = PT.lambda;
                FFsave(a).phimin_Err(l,1) = PT_err(:,1);
                FFsave(a).phimax_Err(l,1) = PT_err(:,2);
                FFsave(a).alpha_Err(l,1) = PT_err(:,4);
                FFsave(a).beta_Err(l,1) = PT_err(:,3);
                FFsave(a).theta_Err(l,1) = PT_err(:,5);
                FFsave(a).lambda_Err(l,1) = PT_err(:,6);
            end
        else
            % Dummies
        end
        if act_chan{a}(5)
            %% Extracting info from RAW structure variable
            for l = 1:numel(target_freq)  % Number of time segments per job
                for w = 2:eigperc  % Eigenvalue percentage
                    Txz(l,w-1) = RAW(w).txz(l,a);
                    Tyz(l,w-1) = RAW(w).tyz(l,a);
                    Txz_Err(l,w-1) = RAW(w).txz_Err(l,a);
                    Tyz_Err(l,w-1) = RAW(w).tyz_Err(l,a);
                end
            end
            %% Polynom fitting
            % txz
            [~,comp,err,~,~,~] = fitpoly2(1./target_freq,real(Txz),eigperc,real(Txz_Err),0,1,polyord);
            FFsave(a).txzr = comp';  % Real
            [~,comp,~,~,~,~] = fitpoly2(1./target_freq,imag(Txz),eigperc,imag(Txz_Err),0,1,polyord);
            FFsave(a).txzi = comp';  % Imaginary
            FFsave(a).txz_Err = err';  % Error
            FFsave(a).txz = complex(FFsave(a).txzr,FFsave(a).txzi);  % Complex
            
            % tyz
            [~,comp,err,~,~,~] = fitpoly2(1./target_freq,real(Tyz),eigperc,real(Tyz_Err),0,1,polyord);
            FFsave(a).tyzr = comp';  % Real
            [~,comp,~,~,~,~] = fitpoly2(1./target_freq,imag(Tyz),eigperc,imag(Tyz_Err),0,1,polyord);
            FFsave(a).tyzi = comp';  % Imaginary
            FFsave(a).tyz_Err = err';  % Error
            FFsave(a).tyz = complex(FFsave(a).tyzr,FFsave(a).tyzi);  % Complex
            
            % Tipper Real and Imaginary Parts
            for i = 1:numel(target_freq)
                FFsave(a).tipr(i) = abs(norm([real(FFsave(a).txz(i)),real(FFsave(a).tyz(i))]));
                FFsave(a).tipi(i) = abs(norm([imag(FFsave(a).txz(i)),imag(FFsave(a).tyz(i))]));
            end
        else
            % Dummies
        end
        FFsave(a).nfreq = numel(target_freq);
        FFsave(a).freq = target_freq';
        FFsave(a).per = 1./(target_freq');
    end
    %% Adding some extra info
    FFsave = FFproc_metainfo(FFsave,job);
       
    %% Saving Eigenvalues information
    EV.best = EV.best(eigperc).best;
    EV.perc = eigperc;
    
    %% Saving results in different folders
    pc = computer;    
    switch pc
        % Windows
        case 'PCWIN64'
            user = getenv('USERPROFILE');
        % Linux    
        case 'GLNXA64'
            user = getenv('USER');
        % MacOS    
        case 'MACI64'
            user = getenv('USER');
    end
    % Path for each FFproc user (Documents)
    cd(fullfile(user,'Documents','FFMT','FFproc','Saving','FFproc Processing'))
    
    % Creating folder name for each job an creating folder (if does not
    % exists)
    file = cell(numel(job.ffts),1);
    for j = 1:numel(job.ffts)
        filename = [job.ffts(j).path,job.ffts(j).file];
        [~,file{j},~] = fileparts(filename);
    end  
    foldername = [];
    if numel(file) == 1
        foldername = file{1};
    else
        for n = 1:numel(job.name)
            foldername = [foldername,[file{n}]];
            if n < numel(file)
                foldername = [foldername,'-'];
            end
        end
    end
    if isfolder(foldername)
    else
        mkdir(foldername)
    end
    cd(foldername)
    % Retrieving target frequency range names
    if target_freq(1) >= 1
        f1 = round(target_freq(1));
        unit1 = 'Hz';
    else
        f1 = round(1./target_freq(1));
        unit1 = 's';
    end
    if target_freq(end) >= 1
        f2 = round(target_freq(end));
        unit2 = 'Hz';
    else
        f2 = round(1./target_freq(end));
        unit2 = 's';
    end
    
    % Outputname
    savename = [num2str(f1),unit1,'-',num2str(f2),unit2,'polyfit_best',num2str(eigperc),'%'];
    
    % Saving depending on variable size
    type = whos('FFsave');
    if type.bytes >= 2*10^9
        save(savename,'FFsave','EV','-v7.3');
    else
        save(savename,'FFsave','EV');
    end

end