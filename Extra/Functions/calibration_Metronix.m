function cal_inst = calibration_Metronix(job,coil,tseg,n)
    
    % Sampling rate
    sr = job.sr;
    calibration = coil_cal_ADU(coil,sr);
    
    % Time interval (Time segment)
    tint = tseg;    
    % Number of samples
    nxt = n;    
    % Frequency vector
    freq = (1:nxt/2)./tint;
    
    cal = cell(1,5);
    for m = 1:3
        coil_amint = interp1(calibration{m}(:,1),calibration{m}(:,2),freq,'pchip');
        coil_phint = interp1(calibration{m}(:,1),calibration{m}(:,3),freq,'pchip');
        calcoil = coil_amint.*exp(1i.*deg2rad(coil_phint));
        cal{m+2} = [calcoil,fliplr(calcoil)];
    end
    
    cal{1} = ones(1,numel(cal{3}));
    cal{2} = cal{1};
    cal_inst = cal;
     
end