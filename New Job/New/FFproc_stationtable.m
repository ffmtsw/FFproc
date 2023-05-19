% function EGstart_stationtable extracts information from ffts global
% variable to be shown in the Station Parameters table.
% version 1.0 / 07abr2020 / cc

function app = FFproc_stationtable(app,ffts)

    if isempty(ffts)
        app.stationtable.Data = [];
    else
        % ID
        ID = num2cell(extractfield(ffts,'ID')');
        % File name
        name = extractfield(ffts,'name')';
        % Device name
        device = extractfield(ffts,'device')';
        % File extension
        ext = extractfield(ffts,'ext')';
        % Sampling rate
        sr = num2cell(extractfield(ffts,'sr')');
        % Start time
        tstart = cellstr(datestr(extractfield(ffts,'tstart')));
        % Stop time
        tend = cellstr(datestr(extractfield(ffts,'tend')));
        % Active channels
        try
            exchan = extractfield(ffts,'exchan')';
            if isempty(exchan)
                exchan = num2cell(logical(ones(numel(ffts),1)));
            end
        catch 
            exchan = num2cell(logical(ones(numel(ffts),1)));
        end
        try
            eychan = extractfield(ffts,'eychan')';
            if isempty(eychan)
                eychan = num2cell(logical(ones(numel(ffts),1)));
            end
        catch 
            eychan = num2cell(logical(ones(numel(ffts),1)));
        end
        try
            bxchan = extractfield(ffts,'bxchan')';
            if isempty(bxchan)
                bxchan = num2cell(logical(ones(numel(ffts),1)));
            end
        catch 
            bxchan = num2cell(logical(ones(numel(ffts),1)));
        end
        try
            bychan = extractfield(ffts,'bychan')';
            if isempty(bychan)
                bychan = num2cell(logical(ones(numel(ffts),1)));
            end
        catch 
            bychan = num2cell(logical(ones(numel(ffts),1)));
        end
        try
            bzchan = extractfield(ffts,'bzchan')';
            if isempty(bzchan)
                bzchan = num2cell(logical(ones(numel(ffts),1)));
            end
        catch 
            bzchan = num2cell(logical(ones(numel(ffts),1)));
        end    
        % Magnetic channels BX
        rsbx = extractfield(ffts,'rsbx')';
        if any(isnan(rsbx)) || isempty(rsbx)
            rsbx = ID;
        else        
            rsbx = num2cell(extractfield(ffts,'rsbx')');
        end
        bxopt = cellstr(string(1:numel(ffts)));
        app.stationtable.ColumnFormat{1,12} = bxopt;    
        % Magnetic channels BY
        rsby = extractfield(ffts,'rsbx')';
        if any(isnan(rsby)) || isempty(rsby)
            rsby = ID;
        else
            rsby = num2cell(extractfield(ffts,'rsby')');         
        end    
        byopt = cellstr(string(1:numel(ffts)));
        app.stationtable.ColumnFormat{1,13} = byopt;
        % Dipoles swap
        try
            exswap = extractfield(ffts,'exswap')';
            if isempty(exswap)
                exswap = num2cell(logical(zeros(numel(ffts),1)));
            end
        catch    
            exswap = num2cell(logical(zeros(numel(ffts),1)));
        end
        try
            eyswap = extractfield(ffts,'eyswap')';
            if isempty(eyswap)
                eyswap = num2cell(logical(zeros(numel(ffts),1)));
            end
        catch    
            eyswap = num2cell(logical(zeros(numel(ffts),1)));
        end
        % Coils swap
        try
            bxswap = extractfield(ffts,'bxswap')';
            if isempty(bxswap)
                bxswap = num2cell(logical(zeros(numel(ffts),1)));
            end
        catch    
            bxswap = num2cell(logical(zeros(numel(ffts),1)));
        end
        try
            byswap = extractfield(ffts,'byswap')';
            if isempty(byswap)
                byswap = num2cell(logical(zeros(numel(ffts),1)));
            end
        catch    
            byswap = num2cell(logical(zeros(numel(ffts),1)));
        end
    
        % Create Table
        fill = [name,...
                device,ext,sr,tstart,tend,...
                exchan,eychan,bxchan,bychan,bzchan,rsbx,rsby,exswap,eyswap,bxswap,byswap];
        
        app.stationtable.Data = fill;
    end
    
end
