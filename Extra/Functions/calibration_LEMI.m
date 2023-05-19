function cal_inst = calibration_LEMI(job,coil,tseg,tlength)
    
    % Sampling rate
    sr = job.sr;   
    % Time interval (Time segment)
    tint = tseg;
    % Number of samples (time segments*sampling rate)
    nxt = tlength;
    % Frequency vector to be evaluated for each target frequency
    freq = (1:nxt/2)./tint;
    
    % Identifying Operative System
    [user,~] = getuser;  
    path = fullfile(user,'Documents','FFMT','FFproc','Calibration Files','LEMI');
    
    if sr > 4 % Calibration for LEMI BBMT systems
        coil = str2double(cellstr(coil));
        % Ex calibration
        cex = ones(1,numel(freq)*2);
        % Ey calibration
        cey = ones(1,numel(freq)*2);
        % Bx calibration
        bxc = load(fullfile(path,[num2str(coil(1)),'.txt']));
        amp = interp1(bxc(:,1),bxc(:,2),freq,'pchip');
        phi = interp1(bxc(:,1),bxc(:,3),freq,'pchip');
        cbx = amp.*exp(1i.*deg2rad(phi));
        cbx = [cbx,fliplr(cbx)];
        % By calibration
        byc = load(fullfile(path,[num2str(coil(2)),'.txt']));
        amp = interp1(byc(:,1),byc(:,2),freq,'pchip');
        phi = interp1(byc(:,1),byc(:,3),freq,'pchip');
        cby = amp.*exp(1i.*deg2rad(phi));
        cby = [cby,fliplr(cby)];
        % Bz calibration
        bzc = load(fullfile(path,[num2str(coil(3)),'.txt']));
        amp = interp1(bzc(:,1),bzc(:,2),freq,'pchip');
        phi = interp1(bzc(:,1),bzc(:,3),freq,'pchip');
        cbz = amp.*exp(1i.*deg2rad(phi));
        cbz = [cbz,fliplr(cbz)];
        
    else % Calibration for LEMI LMT systems (417, 420)
        % Ex calibration
        fid_ex = fopen(fullfile(path,'e1_i.rsp'),'r');
        infoex = textscan(fid_ex,'%f %f %f',inf,'Headerlines',2);
        fclose(fid_ex);
        amp = interp1(infoex{1},infoex{2},freq,'pchip');
        phi = interp1(infoex{1},infoex{3},freq,'pchip');
        cex = amp.*exp(1i.*deg2rad(phi));
        cex = [cex,fliplr(cex)];

        % Ey calibration
        fid_ey = fopen(fullfile(path,'e2_i.rsp'),'r');
        infoey = textscan(fid_ey,'%f %f %f',inf,'Headerlines',2);
        fclose(fid_ey);
        amp = interp1(infoey{1},infoey{2},freq,'pchip');
        phi = interp1(infoey{1},infoey{3},freq,'pchip');
        cey = amp.*exp(1i.*deg2rad(phi));
        cey = [cey,fliplr(cey)];    

        % Bx calibration
        fid_bx = fopen(fullfile(path,'bx_i.rsp'),'r');
        infobx = textscan(fid_bx,'%f %f %f',inf,'Headerlines',2);
        fclose(fid_bx);
        amp = interp1(infobx{1},infobx{2},freq,'pchip');
        phi = interp1(infobx{1},infobx{3},freq,'pchip');
        cbx = amp.*exp(1i.*deg2rad(phi));
        cbx = [cbx,fliplr(cbx)];

        % By calibration
        fid_by = fopen(fullfile(path,'by_i.rsp'),'r');
        infoby = textscan(fid_by,'%f %f %f',inf,'Headerlines',2);
        fclose(fid_by);
        amp = interp1(infoby{1},infoby{2},freq,'pchip');
        phi = interp1(infoby{1},infoby{3},freq,'pchip');
        cby = amp.*exp(1i.*deg2rad(phi));
        cby = [cby,fliplr(cby)];

        % Bz calibration
        fid_bz = fopen(fullfile(path,'bz_i.rsp'),'r');
        infobz = textscan(fid_bz,'%f %f %f',inf,'Headerlines',2);
        fclose(fid_bz);
        amp = interp1(infobz{1},infobz{2},freq,'pchip');
        phi = interp1(infobz{1},infobz{3},freq,'pchip');
        cbz = amp.*exp(1i.*deg2rad(phi));
        cbz = [cbz,fliplr(cbz)];       
    end
        
    % Saving calibration factors
    cal{1} = cex;
    cal{2} = cey;
    cal{3} = cbx;
    cal{4} = cby;
    cal{5} = cbz;
    
    cal_inst = cal;

end