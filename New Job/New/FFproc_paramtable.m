% function FFproc_paramtable extracts information from ffts global variable
% to be shown in the Parameters File. Electrode, Coils, and Coordinate
% tables were concatenated in Parameters table.
% version 1.0 / 07abr2020 / cc
% version 2.0 / 16nov2020 / cc

function FFproc_paramtable(app,ffts)

    if isempty(ffts)
        app.stationtable.Data = [];
    else
        try
            % Name
            name = extractfield(ffts,'name')';
            % Electrodes and dipole rotation
            dip = cat(1,ffts.diplen);
            ex = dip(:,1);
            ey = dip(:,2);
            e_rot = cat(1,ffts.rot);
            % Coil Serial Number and RAP
            coil = cat(1,ffts.coil);
            if isempty(coil)
                coil = {'bx','by','bz'};
                coil = repmat(coil,numel(ffts),1);
            end
            try
                RAP = cat(1,ffts.RAP);
            catch
                RAP = false(numel(ffts),1);
            end            
            % Latitude, longitude and elevation
            coor = cat(1,ffts.GPS);
            lat = coor(:,1);
            lon = coor(:,2);
            z = coor(:,3);
        catch
            % Name
            name = extractfield(ffts,'site')';
            % Electrodes and dipole rotation
            ex = 50*ones(numel(ffts),1);
            ey = 50*ones(numel(ffts),1);
            e_rot = zeros(numel(ffts),1);
            % Coil Serial Number and RAP
            coil = {'bx','by','bz'};
            coil = repmat(coil,numel(ffts),1);
            RAP = false(numel(ffts),1);
            % Latitude, longitude and elevation
            lat = zeros(numel(ffts),1);
            lon = zeros(numel(ffts),1);
            z = zeros(numel(ffts),1);
        end
        
        try
            fill = [name,num2cell(ex),num2cell(ey),num2cell(e_rot),...
                         cellstr(coil),num2cell(RAP),...
                         num2cell(lat),num2cell(lon),num2cell(z)]; 
            app.paramtable.Data = fill;    
        catch
            try
                fill = [name,num2cell(ex),num2cell(ey),num2cell(e_rot),...
                             num2cell(coil),num2cell(RAP),...
                             num2cell(lat),num2cell(lon),num2cell(z)]; 
                app.paramtable.Data = fill;    
            catch
                fill = [name,num2cell(ex),num2cell(ey),num2cell(e_rot),...
                             coil,num2cell(RAP),...
                             num2cell(lat),num2cell(lon),num2cell(z)];
                app.paramtable.Data = fill;    
            end
        end   
    end

end