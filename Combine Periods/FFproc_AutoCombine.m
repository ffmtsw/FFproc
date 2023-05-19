function FFproc_AutoCombine(app)

    path = uigetdir('Select folder which contains all your processed jobs');
    if path == 0
        return
    end

    info = dir(path);

    dlg = uiprogressdlg(app.FFproc,'Title','FFproc','ShowPercentage','on','Cancelable','on');

end