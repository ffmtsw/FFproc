function FFproc_datafromtable(app,ffts)

    station = app.stationtable.Data;
    param = app.paramtable.Data;
    
    % Reasigning data to variables
    for i = 1:numel(ffts)

        % From stationtable
        ffts(i).name = station{i,1};
        ffts(i).tstart = datenum(station{i,5});
        ffts(i).tend = datenum(station{i,6});
        ffts(i).exchan = station{i,7}(1);
        ffts(i).eychan = station{i,8}(1);
        ffts(i).bxchan = station{i,9}(1);
        ffts(i).bychan = station{i,10}(1);
        ffts(i).bzchan = station{i,11}(1);
        ffts(i).rsbx = station{i,12}(1);
        ffts(i).rsby = station{i,13}(1);
        ffts(i).actchan = [ffts(i).exchan,ffts(i).eychan,ffts(i).bxchan,ffts(i).bychan,ffts(i).bzchan];
        ffts(i).exswap = station{i,14}(1);
        ffts(i).eyswap = station{i,15}(1);
        ffts(i).bxswap = station{i,16}(1);
        ffts(i).byswap = station{i,17}(1);

        % From paramtable
        ffts(i).diplen = [param{i,2:3}];
        ffts(i).rot = [param{i,4}];
        ffts(i).coil = param(i,5:7);
        ffts(i).RAP = param{i,8};
        ffts(i).GPS = [param{i,9:end}];
        
    end

    app.ffts = ffts;

end