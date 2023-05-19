function [ffsave,ev] = procsave_load(ffsave,ev)

    [file,path] = uigetfile({'*.mat','FFsave (*.mat)';
                             '(*.*)','All Files (*.*)'},...
                             'Select Process results to combine');    
    cd(path)                 
    % If user press cancel or close button
    if file == 0
        return
    end
    
    load([path,file]);
    if numel(FFsave) > 1
        warning off
        app = procsave_siteselect(FFsave);
        uiwait(app.UIFigure)
        ind = app.ind;
        close(app.UIFigure)
    else
        ind = 1;
    end
            
    FFsave(ind).ID = {file};
    for i = 1:FFsave(ind).nfreq
        FFsave(ind).include(i,1) = true;
        FFsave(ind).marker(i,1) = {'.'};
    end
    
    ev = [ev,EV];
    ffsave = [ffsave,FFsave(ind)];

end
   