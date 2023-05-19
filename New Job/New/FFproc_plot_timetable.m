function FFproc_plot_timetable(app,ffts)

    global font
    if isempty(font)
        font = 'Verdana';
    end
    
    if isempty(ffts)
        cla(app.tsaxes)
        app.tsaxes.Visible = 'off';
    else
        % Time table
        cla(app.tsaxes)
        app.tsaxes.Visible = 'on';
        hold(app.tsaxes,'on')
        for i = 1:length(ffts)
            tbar = ffts(i).t;
            idbar = ones(length(tbar),1).*ffts(i).ID;
            tbar(ffts(i).gaps(1:end-1,2)) = NaN;
            plot(app.tsaxes,datetime(datevec(tbar)),idbar,'linewidth',7)
        end
        grid(app.tsaxes,'on')
    
        if length(ffts) == 1
            set(app.tsaxes,'ytick',1,'yticklabel',1)
        elseif length(ffts) == 2
            set(app.tsaxes,'ytick',1:2,'yticklabel',1:2)
        else
            set(app.tsaxes,'ytick',vertcat(ffts.ID),'yticklabel',vertcat(ffts.ID))
        end
    
        set(app.tsaxes,'YDir','reverse','box','on','Fontsize',10)
        axis(app.tsaxes,'tight')
        xlim(app.tsaxes,datetime(datevec([min(cat(1,ffts.t)) max(cat(1,ffts.t))])))
        ylim(app.tsaxes,[0.5, ffts(i).ID + 0.5])

        xlabel(app.tsaxes,'Date / Time')
        ylabel(app.tsaxes,'Station')
        set(app.tsaxes,'FontName',font,'FontSize',9)
    end

end