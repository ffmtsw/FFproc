% function FFproc_coiltable extracts information from ffts global
% variable to be shown in the Calibration - coils table.
% version 1.0 / 07abr2020 / cc

function app = FFproc_coiltable(app)

    global ffts
    % Coil Serial Number
    coil = (cat(1,ffts.coil));
    RAP = num2cell(cat(1,ffts.RAP));
    if isempty(coil)
        coil = num2cell(zeros(numel(ffts),3));
    end
        
    fill = [coil,RAP];
    set(app.coiltable,'Data',fill)
    
end
