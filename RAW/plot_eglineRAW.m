function plot_eglineRAW(app,info,EV,a,targ)
    
    global font 

    td = tiledlayout(1,1,'Parent',targ,'TileSpacing','compact','Padding','compact');
    ax = nexttile(td);

    cla(ax,'reset')
    hold(ax,'on')
    grid(ax,'on')
    
    per = 1./EV.tfreq;
    best = EV.best(a).best;
 
    plot(ax,per,best(:,1),'color',[1,0.25,0.25],'linewidth',3)
    plot(ax,per,best(:,2),'color',[0.1 0.6 0.9],'linewidth',3)
    plot(ax,per,best(:,3:end),':','color',[.5 .5 .5],'linewidth',0.5)
    
    xlabel(ax,'period (s)')
    ylabel(ax,'EV magnitude')
    title(ax,{'Eigenvalues';['Best-ranked ',num2str(a),'% time segments']})
    legend(ax,'EV1','EV2','other EV','Location','northwest')
    
    set(ax,'Fontname',font,'Fontsize',12,...
           'box','on','Xscale','log','Yscale','log',...
           'xlim',10.^[info.permin info.permax],...
           'xtick',10.^(ceil(info.permin):floor(info.permax)),...
           'MinorGridLineStyle','none')
    
end