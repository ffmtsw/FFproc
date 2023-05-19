function create_cts(app)

    dlg = uiprogressdlg(app.UIFigure,'Title','Creating CTS Files','Indeterminate','on','Cancelable',"on");
    [ren,~] = size(app.table.Data);
    ind_pfc = [];
    for i = 1:ren
        switch app.table.Data{i,6}
            case 'AMT'
                ind_pfc(i,:) = 2:4;
            case 'MT'
                ind_pfc(i,:) = 3:5;
        end
    end

    for j = 1:ren
        % Moving to selected path
        cd(app.table.Data{j,2})
        for i = 1:length(ind_pfc(j,:))
            dlg.Message = ['Creating file: ',app.table.Data{j,1},' - ',app.table.Data{j,i+2}];
            prog_input = ['"',app.syspath.Text,'"'];
            pfc_input = ['"',app.pfcpath.Text,'\PFC_',num2str(ind_pfc(j,i)),'.pfc"'];
            tbl_input = ['"',app.table.Data{j,2},app.table.Data{j,1},'.TBL','"'];
            clb_input = ['"',app.clbpath.Text,'\','*.CLB','"'];
            clc_input = ['"',app.clcpath.Text,'\','*.CLC','"'];
            cts_input = ['"',app.table.Data{j,2},app.table.Data{j,1},'.CTS','"'];
            % Input command for SysCal
            command = [prog_input,' ',pfc_input,' ',tbl_input,' ',clb_input,' ',clc_input,' ',cts_input];
            [status(j,i),result] = system(command);
            % In case of error, continue to the next file
            if status(j,i) == 0 || status(j,i) == 1229
                old = [app.table.Data{j,1},'.CTS'];
                new = [app.table.Data{j,1},'.CTS',num2str(ind_pfc(j,i))];
                movefile(old,new,'f');
            else
            end
        end
        
        % Change status in UITable
        if sum(status(j,:)) == 0 || sum(status(j,:)) == 1229*j
            app.table.Data{j,7} = 'Completed';
        else
            app.table.Data{j,7} = 'Error!';
        end
        
        % Number of succesful jobs
        if isempty(find(status))
            flag(j) = true;
        else
            flag(j) = false;
        end
        
        if dlg.CancelRequested
            uialert(app.figure,'Operation cancelled','CreateCTS - Job Aborted','Icon','Error')
        	return            
        end        
    end
    
    % Output on screen
    beep;
    if sum(flag) == numel(flag)
        icon = 'Success';
    else
        icon = 'Warning';
    end   
    figure(app.UIFigure);
    uialert(app.UIFigure,{['  Correct CTS jobs:',num2str(sum(flag))];...
                          ['  Failed CTS jobs:',num2str(numel(flag)-sum(flag))]},...
                           'CreateCTS - Finish','Icon',icon);

end
