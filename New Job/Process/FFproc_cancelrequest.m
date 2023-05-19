function status = FFproc_cancelrequest(app)

    canreq = uiconfirm(app.UIFigure,'Do you want to cancel the job?','FFproc',...
                                    'Options',{'Current','All','Cancel'},...
                                    'DefaultOption',3,'CancelOption',3,'Icon','warning');

        switch canreq
            case 'Current'
                status = 999;     
            case 'All'
                status = 0;
            case 'Cancel'
                status = NaN;
        end

        app.procprogress.CancelRequested = false;

end