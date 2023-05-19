function FFsave = FFsave2mt(FFsave)
    
    mt = create_mt;
    FFsavenames = fieldnames(FFsave);
    mtnames = fieldnames(mt);    
    comparisonnames = intersect(FFsavenames,mtnames);
    
    for i = 1:numel(FFsavenames)
        if ~ismember(FFsavenames{i},comparisonnames)
                FFsave = rmfield(FFsave,FFsavenames{i});
        end
    end

    comparisonnames = ismember(mtnames,FFsavenames);
    for i = 1:numel(FFsave)
        ffsave = FFsave(i);
        for j = 1:numel(mtnames)
            if comparisonnames(j) == 0
                ffsave = setfield(ffsave,mtnames{j},[]);            
            end            
        end
        MT(i) = ffsave;
        MT(i).ID = i;
        MT(i).source = {'data'};
        MT(i).lonlat = [MT(i).lonlat(2),MT(i).lonlat(1)];
        clear ffsave
    end
    
    FFsave = MT;
    FFsave = orderfields(FFsave,mt); 

end