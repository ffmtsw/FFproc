function resp = read_SR(fr,srfile)

    % Opening calibration file
    fid = fopen(srfile);
    % Scanning data
    cdata = textscan(fid,'%f','delimiter',' ','HeaderLines',2);
    fclose(fid);
    % Find NaNs
    ind = isnan(cdata{1});
    % Info contains only numeric values
    info = cdata{1}(~ind);
    % First 5 values corresponds to factors
    fact = info(1:5);
    % Real and Imaginary coefficients (Ex, Ey, Bx, By, Bz)
    data = reshape(info(6:end),11,[])';
    
    % Frequency vector
    freq = 1./data(:,1);
    % Real part of data
    re = data(:,2:2:end);
    % Imaginary part of data
    im = data(:,3:2:end);  
    
    for i = 1:5
        % Complex response interpolated to target frequencies
        comp = interp1(freq(1:end-1),complex(re((1:end-1),i),im((1:end-1),i)),fr(:));
        resp(:,i) = comp*fact(i);
    end

end