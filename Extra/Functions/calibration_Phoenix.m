function cal_inst = calibration_Phoenix(job,ffts,tseg,n)

    % Sampling rate
    sr = job.sr;
    if sr == 24000
        band = '2';
    elseif sr == 2400
        band = '3';
    elseif sr == 150
        band = '4';
    else
        band = '5';
    end
    
    % Time interval (Time segment)
    tint = tseg;
    % Number of samples 
    nxt = n;
    % Frequency vector
    freq = (1:nxt/2)./tint;
    
    % Creating filename
    [~,filename,~] = fileparts(fullfile(ffts.path,ffts.file));    
    ctsfile = fullfile(ffts.path,[strtok(filename,'_'),'.CTS',band]);
    
    % Extracting coefficients from CTS files
    calibration = read_CTS(freq,ctsfile);
    
    cal = cell(1,5);
    for i = 1:5
        resp = calibration(:,i).';
        cal{i} = [resp,fliplr(resp)];
    end
    cal_inst = cal;

end