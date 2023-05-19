function plot_rhoRAW(app,evperc,RAW,targ)
    
    load('rho1.mat')
    load('rho2.mat')
    
    n = numel(RAW(2).site);
    name = RAW(2).site;
    per = RAW(2).per;
    
    ax = subplot(2,1,1,'Parent',targ);   
    cla(ax)
    hold(ax,'on')
    grid(ax,'on')
    
    for i = 1:n
        scatter(ax,log10(per),log10(RAW(evperc).rhoxy(:,i)),40,ones(numel(per),1)*rho1(i,:),'filled')
        scatter(ax,log10(per),log10(RAW(evperc).rhoyx(:,i)),40,ones(numel(per),1)*rho2(i,:),'filled')
        leg{2*i-1,1} = [name{i},' xy'];
        leg{2*i,1} = [name{i},' yx'];
    end
    xlabel(ax,'log10 period (s)')  
    ylabel(ax,'log10 \rho (\Omegam)')
    title(ax,'Resistivity')

    set(ax,'Fontname','Verdana','Fontsize',12,...
           'box','on')
    
    L = legend(ax,'location','northeast','NumColumns',2,'Fontsize',6);
    L.String = leg;
    

end