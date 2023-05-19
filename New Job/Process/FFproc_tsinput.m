function ts = FFproc_tsinput(path,file,opt)
    
    % Getting information from file
    if numel(string(file)) > 1
        [~,~,ext] = fileparts(fullfile(path,file(1)));  
    else
        [~,~,ext] = fileparts(fullfile(path,file));  
    end
    ext = upper(ext);
    
    readdata = true;
    printdata = false;
    
    % Selection of time series
    switch ext
    % ---------------------------------------------- METRONIX
        case '.ATS'
            ts = read_ATS(path,readdata,printdata,opt);
    % ---------------------------------------------- PHOENIX
        case '.TS2'
            ts = read_TS(path,file,readdata,printdata);
            if isempty(ts)
                return
            end
        case '.TS3'
            ts = read_TS(path,file,readdata,printdata);
            if isempty(ts)
                return
            end            
        case '.TS4'
            ts = read_TS(path,file,readdata,printdata);
            if isempty(ts)
                return
            end
        case '.TS5'
            ts = read_TS(path,file,readdata,printdata);
            if isempty(ts)
                return
            end
    % ---------------------------------------------- LEMI420
        case '.TXT'
            ts = read_LEMI420(path,file,readdata,printdata,opt);
    % ---------------------------------------------- EMERALD
        case '.RAW'
            ts = read_RAW(path,file,readdata,printdata);
    % ---------------------------------------------- GEOLORE
        case '.ML1'
            ts = read_ML1(path,file,readdata,printdata,opt);
        case '.ML2'
            ts = read_ML2(path,file,readdata,printdata,opt);
    % --------------------------------------------- Intermagnet Observatory
        case '.SEC'
            ts = read_SEC(path,file,readdata,printdata);
        case '.MIN'
            ts = read_SEC(path,file,readdata,printdata);
    % ---------------------------------------------- MatlabTS       
        case '.MAT'
            try 
                % Exported from ReadTS
                ts = read_MAT(path,file,readdata,printdata);
            catch
                try
                    % CR Device (Nawa, 2020)
                    ts = read_CR(path,file,readdata,printdata);
                catch
                    % Artificial Time Series
                    ts = read_ART(path,file,readdata,printdata);
                end
            end            
    % ------------------------------------------------------------- UPPSALA
        case '.ASC'
            ts  = read_ASC(path,file,readdata,printdata);            
            
        otherwise
            aux = ext;
            ext = regexprep(ext,'[\d"]','');            
            switch ext
             % ---------------------------------------------- LEMI
                case '.T'
                    ts = read_LEMI417(path,file,readdata,false,opt);
                case '.B'
                    ts = read_B(path,file,readdata,false,opt);
            end
    end
    
end