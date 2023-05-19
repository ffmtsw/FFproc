function plot_egbars(app,EV,scale,perc)

    tfreq = extractfield(EV,'tfreq')';
    ntseg = extractfield(EV,'ntseg')';
    EVal = extractfield(EV,'EVal')';
    EVal = horzcat(EVal{:});
    timeseg = extractfield(EV,'timeseg')';
    timeseg = horzcat(timeseg{:});
    barw = 0.02;
    
    delete(findall(app.evbargrid,'type','uitabgroup'))
    app.tabgroup = uitabgroup(app.evbargrid,'TabLocation','bottom');    

    for i = 1:length(tfreq)        
        % EV        
        ev1 = abs(EVal(i).EVal(:,1));
        ev2 = abs(EVal(i).EVal(:,2));
        ev3 = abs(EVal(i).EVal(:,3));      
        % Bins EV1
        bmin = floor(min(log10(ev1))/barw)*barw;
        bmax = ceil(max(log10(ev1))/barw)*barw;
        bev1 = bmin:barw:bmax;
        % Bins EV2
        bmin = floor(min(log10(ev2))/barw)*barw;
        bmax = ceil(max(log10(ev2))/barw)*barw;
        bev2 = bmin:barw:bmax;
        % Bins EV3
        bmin = floor(min(log10(ev3))/barw)*barw;
        bmax = ceil(max(log10(ev3))/barw)*barw;
        bev3 = bmin:barw:bmax;
        
        % EV after threshold selection
        ev1_sel = ev1(timeseg(i).timeseg(1:ceil(perc*ntseg(i)/100)));
        ev2_sel = ev2(timeseg(i).timeseg(1:ceil(perc*ntseg(i)/100)));
        ev3_sel = ev3(timeseg(i).timeseg(1:ceil(perc*ntseg(i)/100))); 
        % Bins EV1 after threshold selection
        bmin = floor(min(log10(ev1_sel))/barw)*barw;
        bmax = ceil(max(log10(ev1_sel))/barw)*barw;
        bev1_sel = bmin:barw:bmax;
        % Bins EV2 after threshold selection
        bmin = floor(min(log10(ev2_sel))/barw)*barw;
        bmax = ceil(max(log10(ev2_sel))/barw)*barw;
        bev2_sel = bmin:barw:bmax;
        % Bins EV3 after threshold selection
        bmin = floor(min(log10(ev3_sel))/barw)*barw;
        bmax = ceil(max(log10(ev3_sel))/barw)*barw;
        bev3_sel = bmin:barw:bmax;
        
        % Histogram calculation
        evbar1 = hist(log10(ev1),bev1);
        evbar2 = hist(log10(ev2),bev2);
        evbar3 = hist(log10(ev3),bev3);
        evbar1_sel = hist(log10(ev1_sel),bev1_sel);
        evbar2_sel = hist(log10(ev2_sel),bev2_sel);
        evbar3_sel = hist(log10(ev3_sel),bev3_sel);
        
        % Plots
        tab(i) = uitab('Parent',app.tabgroup,'BackgroundColor',[1 1 1],...
                       'Tag',num2str(i),'AutoResizeChildren','off');
        if tfreq(i) > 1
            tab(i).Title = [num2str(round(tfreq(i))),' [Hz]'];
        else
            tab(i).Title = [num2str(round(1./tfreq(i))),' [s]'];
        end        
        td = tiledlayout(1,3,'Parent',tab(i),'TileSpacing','compact','Padding','compact');

        % Eigenvalue 1
        ax1 = nexttile(td);
        app.evbar1(i) = ax1;
        barh(ax1,bev1,evbar1,'Edgecolor','none','Facecolor',[0.8,0.8,0.8]);   
        hold(ax1,'on')
        barh(ax1,bev1_sel,evbar1_sel,'Edgecolor','none','Facecolor',[0.47,0.67,0.19]);
        if bev3(1) > 0
            ylim(ax1,[0.5*bev3(1) 1.1*bev1(end)])
        else
            ylim(ax1,[1.5*bev3(1) 1.1*bev1(end)])
        end
        xlim(ax1,[0 1.2*max(evbar1)])
        title(ax1,'Eigenvalue 1')
        xlabel(ax1,'counts')
        ylabel(ax1,'log10 signal/noise')
        tick = get(ax1,'Yticklabel');
        tick = round(10.^(str2double(tick))/0.01)*0.01;
        set(ax1,'XScale',scale,'Fontname','Verdana','Fontsize',14)
        grid(ax1,'on')
        
        % Eigenvalue 2
        ax2 = nexttile(td);
        app.evbar2(i) = ax2;
        barh(ax2,bev2,evbar2,'Edgecolor','none','Facecolor',[0.8,0.8,0.8]);              
        hold(ax2,'on')
        barh(ax2,bev2_sel,evbar2_sel,'Edgecolor','none','Facecolor',[0.47,0.67,0.19]);
        if bev3(1) > 0
            ylim(ax2,[0.5*bev3(1) 1.1*bev1(end)])
        else
            ylim(ax2,[1.5*bev3(1) 1.1*bev1(end)])
        end
        xlim(ax2,[0 1.2*max(evbar2)])
        title(ax2,'Eigenvalue 2')
        xlabel(ax2,'counts')
        tick = get(ax2,'Yticklabel');
        tick = round(10.^(str2double(tick))/0.01)*0.01;
        set(ax2,'XScale',scale,'Fontname','Verdana','Fontsize',14)
        grid(ax2,'on')
        
        % Eigenvalue 3
        ax3 = nexttile(td);
        app.evbar3(i) = ax3;
        barh(ax3,bev3,evbar3,'Edgecolor','none','Facecolor',[0.8,0.8,0.8]); 
        hold(ax3,'on')
        barh(ax3,bev3_sel,evbar3_sel,'Edgecolor','none','Facecolor',[0.47,0.67,0.19]);
        if bev3(1) > 0
            ylim(ax3,[0.5*bev3(1) 1.1*bev1(end)])
        else
            ylim(ax3,[1.5*bev3(1) 1.1*bev1(end)])
        end
        xlim(ax3,[0 1.2*max(evbar3)])
        title(ax3,'Eigenvalue 3')
        xlabel(ax3,'counts')
        tick = get(ax3,'Yticklabel');
        tick = round(10.^(str2double(tick))/0.01)*0.01;
        set(ax3,'XScale',scale,'Fontname','Verdana','Fontsize',14)
        grid(ax3,'on')        

        lg = legend(ax3,'All data',['EV criteria selected ',num2str(perc),'%'],'Fontsize',10,'Orientation','Vertical','NumColumns',2);
        lg.Layout.Tile = 'South';   
    end

    app.tab = tab;

end