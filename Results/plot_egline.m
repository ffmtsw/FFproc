function plot_egline(info,FFsave,EV,targ,pos)
    
    global font 
    if isempty(pos)
        ax = nexttile(targ);
    else
        ax = nexttile(targ,pos); 
    end  
    cla(ax,'reset')
    hold(ax,'on')
    grid(ax,'on')
    
    per = FFsave.per;
    best = cat(1,EV.best);    
    [per,id] = sort(per);
    
    plot(ax,per,best(id,1),'color',[1,0.25,0.25],'linewidth',3)
    plot(ax,per,best(id,2),'color',[0.1 0.6 0.9],'linewidth',3)
    plot(ax,per,best(id,3:end),':','color',[.5 .5 .5],'linewidth',0.5) 
    
    xlabel(ax,'period (s)')
    ylabel(ax,'EV magnitude')
    title(ax,['Best',num2str(EV(end).perc),'% Ranked Eigenvalues']) 
    legend(ax,'EV1','EV2','other EV','Location','northwest')
    
    set(ax,'Fontname',font,'Fontsize',12,...
           'box','on','Xscale','log','Yscale','log',...
           'xlim',10.^[info.permin info.permax],...
           'xtick',10.^(ceil(info.permin):floor(info.permax)),...
           'MinorGridLineStyle','none')
       
   app.eglineax = ax;
    
end