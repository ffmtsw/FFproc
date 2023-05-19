function RAW = FFproc_RAW_metainfo(RAW,job)

    for i = 1:numel(RAW)
        for j = 1:numel(job.ffts)
            
            % Coordinates
            RAW(i).site{j} = job.name{j};
            [x,y,~] = deg2utm(job.ffts(j).GPS(1),job.ffts(j).GPS(2));
            RAW(i).lonlat(j,:) = [job.ffts(j).GPS(1),job.ffts(j).GPS(2)];
            RAW(i).UTM(j,:) = [x,y];
            RAW(i).z(j) = job.ffts(j).GPS(3);
            
        end
    end

end