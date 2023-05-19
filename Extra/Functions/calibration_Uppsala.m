function cal_inst = calibration_Uppsala(job,ffts,tseg)

    sr = job.sr;
    tint = tseg;
    nxt = floor(tint*sr);
    
    freq = (1:nxt/2)./tint;
    
    % Identifying Operative System
    [user,~] = getuser;
    
    % Name of calibration file
    name = 'met_coils_010203.sr';
    srfile = fullfile(user,'Documents','FFMT','FFproc','Calibration Files','Uppsala',name);
    
    % Retrieving coeficients from calibration files
    calibration = read_SR(freq,srfile);

    cal = cell(1,5);
    for i = 1:5
        resp = calibration(:,i).';
        cal{i} = [resp,fliplr(resp)];
    end
    cal_inst = cal;

end