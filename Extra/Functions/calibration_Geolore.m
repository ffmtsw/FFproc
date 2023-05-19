function cal_inst = calibration_Geolore(job,coil,tseg,n,RAP)          
   
    % Time interval (Time segment)
    tint = tseg;
    % Number of samples 
    nxt = n;
    % Frequency vector
    freq = (1:nxt/2)./tint;
    % Angular frequency
    omegacal = 2*pi*freq;
    
    % Identifying Operative System
    [user,~] = getuser;
    
    % Geolore Electronic components
    R1 = 47*1e3;
    R3 = 22*1e3;
    C1 = 2.2/1e6;
    C3 = 4.7/1e6; 
    
    % Active filter
    F1 = 1./(1+(1i*C1.*omegacal).*(2*R1)-omegacal.^2.*(R1*R1*C1*C1));
    % Passive filter
    F2 = 1./(1+1i*omegacal*R3*C3); 
    
    % Passive one pole filter (Only for Geolore 3rd Generation)
%     if job.sr == 10
    if strcmp(coil{1}(1:2),'AU')
        F3 = 1./(1 + (1i*omegacal*100*1000*(0.1/1000000)));
    else
        F3 = 1;
    end
    % Electric field filter
    Fmult = F1.*F2.*F3;
    
    F_el = [Fmult,fliplr(Fmult)];
    
    % Geomag calibration
    if (isnumeric(str2double(coil{1})) || isnumeric(str2double(coil(1)))) ...
        && ~isnan(str2double(coil{1}))      
        fid = fopen(fullfile(user,'Documents','FFMT','FFproc', ...
            'Calibration Files','GEOLORE','MagFil.txt'),'r');
        if fid~=(-1)
            mfil = (fscanf(fid,'%f',[20 inf]));
        fclose(fid);
        end
        amint = interp1(mfil(:,1),mfil(:,2),omegacal./(2*pi),'pchip');
        phint = interp1(mfil(:,1),mfil(:,3),omegacal./(2*pi),'pchip');

        % Magnetic field filter
        F_mag = amint.*exp(1i.*(phint.*(pi/180)));
        F_mag = 0.75*[F_mag,fliplr(F_mag)];

        cal{1} = F_el;
        cal{2} = F_el;
        cal{3} = F_mag;
        cal{4} = F_mag;
        cal{5} = F_mag;
            
    % Magson Calibration
    elseif strcmp(coil{1}(1:2),'AU')        
        k = [0.5847 1.4706 0.2711;...
             0.8033 0.9494 0.3722;...
             1.1756 0.5130 0.1764;...
             0.0000 0.9449 0.0000];

        for i = 1:4
            A(:,i) = 1-((omegacal.^2).*k(i,1).*k(i,2));
            B(:,i) = omegacal.*(k(i,2)+k(i,3));
            F(:,i) = 1./(A(:,i) + 1i*B(:,i));
        end
        
        F_mag = prod(F,2);
%         F_mag = 0.75*[F_mag',fliplr(F_mag')];
        F_mag = [F_mag',fliplr(F_mag')];
        
        if RAP
            load(fullfile(user,'Documents','FFMT','FFproc','Calibration Files','Geolore','MAGSON_RAP_with.mat'));
        else
            load(fullfile(user,'Documents','FFMT','FFproc','Calibration Files','Geolore','MAGSON_RAP_without.mat'));
        end
        
        ind = find(strcmp(coil(1),Coef.MAG));
        coefx = Coef.X(ind);
        coefy = Coef.Y(ind);
        coefz = Coef.Z(ind); 
        
        cal{1} = 1./F_el;
        cal{2} = 1./F_el;
        cal{3} = 1./(F_mag.*coefx);
        cal{4} = 1./(F_mag.*coefy);
        cal{5} = 1./(F_mag.*coefz);
    end
        
    cal_inst = cal;

end