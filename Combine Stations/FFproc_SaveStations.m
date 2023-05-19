function FFproc_SaveStations(app)
    
    % Retrieving information from input stations
    for i = 1:size(app.outtable.Data,1)
        info(i).folder = app.outtable.Data{i,1};
        info(i).name = app.outtable.Data{i,2};
    end
    
    % Prealocating space for ffsave and ev variables
    MT = [];
%     ev = [];
    for i = 1:numel(info)
        load(fullfile(info(i).folder,info(i).name))
        if exist("FFsave",'var')
            mt = FFsave;
        end
        MT = [mt,MT];        
%         if i > 1
%             if isequal(ev,EV)
%                 % Updating EV variable
%                 ev = EV;
%             else
%                 answer = uiconfirm(app.figure,{'EV structures do not match';...
%                                                'EV structure will be empty';...
%                                                'Continue?'},...
%                                                'FFproc - Combine Stations',...
%                                                'Options',{'Yes','No'},...
%                                                'DefaultOption',1,'CancelOption',2,...
%                                                'Icon','Warning');
%                 switch answer
%                     case 'Yes'
%                         % Removing ev variable
%                         ev = [];
% %                         continue
%                     case 'No'
%                         uialert(app.figure,'Remove incompatible site(s)',...
%                                            'FFproc - Saving','Icon','Warning');
%                         return
%                 end
%             end
%         else
%             % Updating EV variable
%             ev = EV;
%         end     
        
    end

    % Updating values for FFsave and EV variables
    mt = MT;
    for i = 1:numel(mt)
        mt(i).ID = i;
    end
%     EV = ev;
   
    %% Saving file
    [user,~] = getuser;
    % Path for each FFproc user (Documents)
    cd(fullfile(user,'Documents','FFMT','FFproc','Saving','FFproc2MTDIM'))    
    
    % Saving data in mat file
    [file,path] = uiputfile('*.mat','Save combined stations','.mat');
    
    try
        type = whos('mt');
        if type.bytes >= 2*10^9
            save(fullfile(path,file),'mt','-v7.3');
        else
            save(fullfile(path,file),'mt');
        end

        uialert(app.UIFigure,{['File: ',file];'saved in:';path},...
                             'Saving','Icon','Success');
    catch
        uialert(app.UIFigure,{['File: ',file];'was not saved!'},...
                              'Saving','Icon','Error');
    end

end