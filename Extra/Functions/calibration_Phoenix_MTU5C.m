function cal = calibration_Phoenix_MTU5C(job,ffts,tseg,n)
   
    % Frequency vector
    freq = (1:n/2)./tseg;
    
    % Extract file names
    info = dir('*.json');
    
    % Extract info from config and recmeta files
    json = readjson(ffts.path);
    
    %% ----------------------- Receiver responses ----------------------- %
    chans = json.recmeta.chconfig.chans;
    
    % Extracting information from channel configuration
    ch = {'E1','E2','H1','H2','H3'};
    % Pre-allocate memory
    ga = NaN(1,numel(chans));
    pg = NaN(1,numel(chans));
    lp = NaN(1,numel(chans));
    at = NaN(1,numel(chans));
    for i = 1:numel(ch)
        for j = 1:numel(chans)
            if strcmp(ch{i},chans{j}.tag)
                ga(i) = chans{j}.ga;    % Gain
                pg(i) = chans{j}.pg;    % Gain preamplifier
                lp(i) = chans{j}.lp;    % Low-pass filter
                at(i) = chans{j}.at;    % Attenuator
            end
        end
    end
    
    % Selecting for correct receiver response based on Gain and LPF
    % Pre-allocate memory
    rx_id = NaN(1,numel(chans));
    switch strtok(json.recmeta.chconfig.chans{1}.hw,'-')
        case 'BCM01' % Before 2022
            for i = 1:numel(ch)
                if at(i) == 0
                    Gain = pg(i)*ga(i);
                    LPF = lp(i);
                else
                    Gain = pg(i)*ga(i)*0.1;
                    LPF = lp(i);
                end      
                % Find the correct response
                if Gain == 1 && LPF == 10000
                    rx_id(i) = 1;
                elseif Gain == 4 && LPF == 10000
                    rx_id(i) = 2;
                elseif Gain == 16 && LPF == 10000
                    rx_id(i) = 3;
                elseif Gain == 1 && LPF == 1000
                    rx_id(i) = 4;
                end
            end
        case 'BCM05' % From 2022
            for i = 1:numel(ch)
                if i == 1 || i == 2     % E channels
                    pg(i) = 8;
                else                    % B Channels
                    pg(i) = 1;
                end
                if at(i) == 0                    
                    Gain = pg(i)*ga(i);
                    LPF = lp(i);
                else
                    Gain = pg(i)*ga(i)*0.1;
                    LPF = lp(i);
                end      
                % Find the correct response
                if Gain == 4 && LPF == 10000
                    rx_id(i) = 1;
                elseif Gain == 8 && LPF == 10000
                    rx_id(i) = 2;
                elseif Gain == 16 && LPF == 10000
                    rx_id(i) = 3;
                elseif Gain == 32 && LPF == 1000
                    rx_id(i) = 4;
                end   
            end
    end
    
    % Opening receiver response file in json format
    rsxid = find(str2double(json.recmeta.instid) == str2double(strtok(extractfield(info,'name'),'_')));
    rsxcal = extractjson(fullfile(ffts.path,info(rsxid).name));
    for i = 1:numel(ch)
        ch_id = find(strcmp(ch(i),extractfield(rsxcal.cal_data,'tag')));
        rsx(i).f = rsxcal.cal_data(ch_id).chan_data(rx_id(i)).freq_Hz;
        rsx(i).amp = rsxcal.cal_data(ch_id).chan_data(rx_id(i)).magnitude;
        rsx(i).phi = rsxcal.cal_data(ch_id).chan_data(rx_id(i)).phs_deg;
    end
    
    % EX receiver response correction
    aresp = interp1(rsx(1).f,rsx(1).amp,freq,'pchip');
    presp = interp1(rsx(1).f,rsx(1).phi,freq,'pchip');
    ex_r = aresp.*exp(1i.*deg2rad(presp));
    % EY receiver response correction
    aresp = interp1(rsx(2).f,rsx(2).amp,freq,'pchip');
    presp = interp1(rsx(2).f,rsx(2).phi,freq,'pchip');
    ey_r = aresp.*exp(1i.*deg2rad(presp));
    % BX receiver response correction
    aresp = interp1(rsx(3).f,rsx(3).amp,freq,'pchip');
    presp = interp1(rsx(3).f,rsx(3).phi,freq,'pchip');
    bx_r = aresp.*exp(1i.*deg2rad(presp));
    % BY receiver response correction
    aresp = interp1(rsx(4).f,rsx(4).amp,freq,'pchip');
    presp = interp1(rsx(4).f,rsx(4).phi,freq,'pchip');
    by_r = aresp.*exp(1i.*deg2rad(presp));
    % BZ receiver response correction
    aresp = interp1(rsx(5).f,rsx(5).amp,freq,'pchip');
    presp = interp1(rsx(5).f,rsx(5).phi,freq,'pchip');
    bz_r = aresp.*exp(1i.*deg2rad(presp));
    
    %% ------------------------ Sensor responses ------------------------ %    
    % BX Magnetic sensor response correction 
    coil = find(ffts.coil{1} == str2double(strtok(extractfield(info,'name'),'_')));  
    BX = extractjson(fullfile(ffts.path,info(coil).name));
    % Interpolating to target frequency
    f = BX.cal_data.chan_data.freq_Hz;
    amp = BX.cal_data.chan_data.magnitude;
    phi = BX.cal_data.chan_data.phs_deg;
    aresp = interp1(f,amp,freq,'pchip');
    presp = interp1(f,phi,freq,'pchip');
    bx_s = aresp.*exp(1i.*deg2rad(presp));

    % BY Magnetic sensor response correction 
    coil = find(ffts.coil{2} == str2double(strtok(extractfield(info,'name'),'_')));    
    BY = extractjson(fullfile(ffts.path,info(coil).name));
    % Interpolating to target frequency
    f = BY.cal_data.chan_data.freq_Hz;
    amp = BY.cal_data.chan_data.magnitude;
    phi = BY.cal_data.chan_data.phs_deg;
    aresp = interp1(f,amp,freq,'pchip');
    presp = interp1(f,phi,freq,'pchip');
    by_s = aresp.*exp(1i.*deg2rad(presp));

    % BZ Magnetic sensor response correction 
    coil = find(ffts.coil{3} == str2double(strtok(extractfield(info,'name'),'_')));    
    BZ = extractjson(fullfile(ffts.path,info(coil).name));
    % Interpolating to target frequency
    f = BZ.cal_data.chan_data.freq_Hz;
    amp = BZ.cal_data.chan_data.magnitude;
    phi = BZ.cal_data.chan_data.phs_deg;
    aresp = interp1(f,amp,freq,'pchip');
    presp = interp1(f,phi,freq,'pchip');
    bz_s = aresp.*exp(1i.*deg2rad(presp));

    %% Final coefficients
    cal = cell(1,5);

    cal{1} = [ex_r;flipud(ex_r)];
    cal{2} = [ey_r;flipud(ey_r)];
    cal{3} = [bx_r.*bx_s;flipud(bx_r.*bx_s)];
    cal{4} = [by_r.*by_s;flipud(by_r.*by_s)];
    cal{5} = [bz_r.*bz_s;flipud(bz_r.*bz_s)];

end