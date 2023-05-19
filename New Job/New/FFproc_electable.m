% function FFproc_electable extracts information from ffts global
% variable to be shown in the Calibration - electrodes table.
% version 1.0 / 07abr2020 / cc

function app = FFproc_electable(app)

    global ffts
    % Electrodes and dipole rotation
    try
        dip = cat(1,ffts.diplen);
        ex_l = num2cell(dip(:,1));
        ey_l = num2cell(dip(:,2));
        e_rot = num2cell(cat(1,ffts.rot));
    catch
        ex_l = num2cell(50*ones(numel(ffts),1));
        ey_l = num2cell(50*ones(numel(ffts),1));
        e_rot = num2cell(zeros(numel(ffts),1));
    end
    fill = [ex_l,ey_l,e_rot];
    
    set(app.electable,'Data',fill)
    
end
