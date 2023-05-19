function procsave_filltable(app,mt)
     
    vect = [];
    for i = 1:numel(mt)
        vect = [vect;ones(numel(mt(i).per),1)*i];
    end
    vect = cellstr(num2str(vect));
    
    per = cat(1,mt.per);

    % Selecting the right label (scale Hz or s)
    for i = 1:numel(per)
        if per(i) < 1
            mag = '(Hz)';
            pernum = 1/per(i);
        else
            mag = '(s)';
            pernum = per(i);
        end
        perlist{i,1} = [num2str(round(pernum,4)),' ',mag];
    end

    include = num2cell(logical(ones(numel(vect),1)));
    fill = [vect,perlist,include];
        
    app.table.Data = fill;
    
end