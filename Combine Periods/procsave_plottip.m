function procsave_plottip(app,mt,info)
   
    global font
    ind = find(cat(1,app.table.Data{:,3}));
    per = cat(1,mt.per);
    txz = cat(1,mt.txz);
    tyz = cat(1,mt.tyz);
    
    ax = subplot(1,3,1,'Parent',app.Panel);
    cla(ax)
    hold(ax,'on')
    grid(ax,'on')
    
    quiver(ax,0,log10(per(ind)),real(tyz(ind)),-real(txz(ind)),...
             'color','r','linewidth',1,'maxheadsize',0.1,'autoscale','off')
    
    quiver(ax,0,log10(per(ind)),imag(tyz(ind)),-imag(txz(ind)),...
             'color','b','linewidth',1,'maxheadsize',0.1,'autoscale','off')
       
    set(ax,'Fontname',font,'Fontsize',10,...
           'box','on',...
           'ydir','reverse',...
           'dataaspectratio',[1 1 1],...
           'xlim',[-0.6 0.6],'xtick',-5:0.5:5,...
           'ylim',[info.permin info.permax])
       
    xlabel(ax,'amplitude')  
    ylabel(ax,'log10 period (s)')
    title(ax,'Tipper')
    
end