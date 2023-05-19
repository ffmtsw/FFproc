function FFsave = FFproc_OPT_metainfo(FFsave,job)

    for i = 1:numel(FFsave)     
        
        FFsave(i).site = job.name(i);
        % Coordinates
        try
            lat = job.ffts(i).GPS(2);
            lon = job.ffts(i).GPS(1);
            z = job.ffts(i).GPS(3); 
            [x,y,~] = deg2utm(lat,lon);
        catch
            lat = job.ffts(i).GPS{1}(i,2);
            lon = job.ffts(i).GPS{1}(i,1);
            z = job.ffts(i).GPS{2}(i);
            [x,y,~] = deg2utm(lat,lon);
        end        
        FFsave(i).lonlat = [lon,lat];
        FFsave(i).UTM = [x,y];
        FFsave(i).z = z; 
        
    end

end