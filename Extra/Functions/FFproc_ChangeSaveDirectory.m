function FFproc_ChangeSaveDirectory(app)

    load('ffprocsavepath.mat','path')
    answer = uiconfirm(app.FFproc,{'Current directory:';'';path;'';'Do you want to change it?'},...
                                  'FFproc','Options',{'Continue','Default','Cancel'},...
                                  'DefaultOption',1,'CancelOption',3,'Icon','info');
    switch answer
        case 'Continue'
            pathnew = uigetdir(pwd,'Select directory to  save FFproc processing results');
            if pathnew == 0
                return
            end
            path = pathnew;
            [user,~] = getuser;
            pathsave = fullfile(user,'Documents','FFMT','Extras');
            save(fullfile(pathsave,'ffprocsavepath.mat'),'path')
            uialert(app.FFproc,{'New saving directory:';'';path;''},'FFproc','Icon','success')

        case 'Default'
            [user,~] = getuser;
            % Path for each FFproc user (Documents)
            path = fullfile(user,'Documents','FFMT','FFproc','Saving','FFproc Processing');
            pathsave = fullfile(user,'Documents','FFMT','Extras');
            save(fullfile(pathsave,'ffprocsavepath.mat'),'path')
            uialert(app.FFproc,{'New saving directory:';'';path;''},'FFproc','Icon','success')
            
        otherwise
            
            return
    end

end