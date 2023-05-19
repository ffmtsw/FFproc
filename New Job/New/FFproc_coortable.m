% function FFproc_coortable extracts information from ffts global
% variable to be shown in the Calibration - coordinates table.
% version 1.0 / 07abr2020 / cc

function app = FFproc_coortable(app)

    global ffts
    % Latitude, longitude and elevation
    try
        coor = cat(1,ffts.GPS);
        lat = num2cell(coor(:,1));
        lon = num2cell(coor(:,2));
        z = num2cell(coor(:,3));
    catch
        lat = num2cell(zeros(numel(ffts),1));
        lon = num2cell(zeros(numel(ffts),1));
        z = num2cell(zeros(numel(ffts),1));
    end
    
    fill = [lat,lon,z];
    
    set(app.coortable,'Data',fill)
    
end
