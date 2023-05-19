function plot_Zev(app,info,RAW,ranking,targ,id,cf)
    
    warning off   
    per = RAW(1).per;
    name = RAW(1).site;

    % Retrieving all estimations
    pertot = cat(1,RAW(1:ranking).per);

    ZXY = cat(1,RAW(1:ranking).Zxy);
    ZXY = reshape(ZXY(:,app.a),[],ranking);

    ZYX = cat(1,RAW(1:ranking).Zyx);
    ZYX = reshape(ZYX(:,app.a),[],ranking);

    ind = repelem((1:ranking)',numel(per),1);
    
    % Extending the period vector
    delta1 = abs(per(2)-per(1));
    delta2 = abs(per(end)-per(end-1));
    perdelta = [per(1)+10^log10(delta1);per;per(end)-10^log10(delta2)];

    td = tiledlayout(2,2,'Parent',targ,'TileSpacing','compact','Padding','compact');
    
    % ------------------------------------------------------------ Real Zxy
    ax = nexttile(td);   
    cla(ax,'reset')
    hold(ax,'on')
    grid(ax,'on')       
    colormap(ax,colors('viridis',ranking-1))   

    scatter(ax,pertot,log10(real(ZXY(:))),50,ind,'filled')
    if exist('id','var') && ~isempty(id)        
        plot(ax,perdelta,cf(app.a).Zxy_Re(log10(perdelta)),'k:','LineWidth',1.5)
        sel = id(app.a).Zxy_Re(:);
        ZXY_sel = NaN(size(sel));
        for i = 1:numel(sel)
            ZXY_sel(i) = ZXY(i,sel(i));
        end
        scatter(ax,per,log10(real(ZXY_sel)),60,sel,'filled','MarkerEdgeColor','k','LineWidth',1.5)
    end 
    xlabel(ax,'period (s)')  
    ylabel(ax,'\Re log10 Zxy (mV/Km/nT)')
    title(ax,['\Re Zxy ',name(app.a)])   
    
    set(ax,'Fontname','Verdana','Fontsize',12,...
           'box','on','Xscale','log','Yscale','linear',...
           'xlim',10.^[info.permin info.permax],...
           'xtick',10.^(ceil(info.permin):floor(info.permax)),...
           'MinorGridLineStyle','none')
    
    % ------------------------------------------------------- Imaginary Zxy
    ax = nexttile(td);     
    cla(ax,'reset')
    hold(ax,'on')
    grid(ax,'on')        
    colormap(ax,colors('viridis',ranking-1)) 
    
    scatter(ax,pertot,log10(imag(ZXY(:))),50,ind,'filled')
    if exist('id','var') && ~isempty(id)        
        plot(ax,perdelta,cf(app.a).Zxy_Im(log10(perdelta)),'k:','LineWidth',1.5)
        sel = id(app.a).Zxy_Im(:);
        ZXY_sel = NaN(size(sel));
        for i = 1:numel(sel)
            ZXY_sel(i) = ZXY(i,sel(i));
        end
        scatter(ax,per,log10(imag(ZXY_sel)),60,sel,'filled','MarkerEdgeColor','k','LineWidth',1.5)
    end 
    xlabel(ax,'period (s)')  
    ylabel(ax,'\Im log10 Zxy (mV/Km/nT)')
    title(ax,['\Im Zxy ',name(app.a)])    

    set(ax,'Fontname','Verdana','Fontsize',12,...
           'box','on','Xscale','log','Yscale','linear',...
           'xlim',10.^[info.permin info.permax],...
           'xtick',10.^(ceil(info.permin):floor(info.permax)),...
           'MinorGridLineStyle','none')
    
    % ------------------------------------------------------------ Real Zyx
    ax = nexttile(td);   
    cla(ax,'reset')
    hold(ax,'on')
    grid(ax,'on')    
    colormap(ax,colors('viridis',ranking-1)) 
    
    scatter(ax,pertot,log10(real(ZYX(:))),50,ind,'filled')
    if exist('id','var') && ~isempty(id)        
        plot(ax,perdelta,cf(app.a).Zyx_Re(log10(perdelta)),'k:','LineWidth',1.5)
        sel = id(app.a).Zyx_Re(:);
        ZYX_sel = NaN(size(sel));
        for i = 1:numel(sel)
            ZYX_sel(i) = ZYX(i,sel(i));
        end
        scatter(ax,per,log10(real(ZYX_sel)),60,sel,'filled','MarkerEdgeColor','k','LineWidth',1.5)
    end 
    xlabel(ax,'period (s)')
    ylabel(ax,'\Re log10 Zyx (mV/Km/nT)')
    title(ax,['\Re Zyx ',name(app.a)])        

    set(ax,'Fontname','Verdana','Fontsize',12,...
           'box','on','Xscale','log','Yscale','linear',...
           'xlim',10.^[info.permin info.permax],...
           'xtick',10.^(ceil(info.permin):floor(info.permax)),...
           'MinorGridLineStyle','none')
    
    % ------------------------------------------------------- Imaginary Zyx
    ax = nexttile(td);   
    cla(ax,'reset')
    hold(ax,'on')
    grid(ax,'on')    
    colormap(ax,colors('viridis',ranking-1)) 

    scatter(ax,pertot,log10(imag(ZYX(:))),50,ind,'filled')
    if exist('id','var') && ~isempty(id)        
        plot(ax,perdelta,cf(app.a).Zyx_Im(log10(perdelta)),'k:','LineWidth',1.5)
        sel = id(app.a).Zyx_Im(:);
        ZYX_sel = NaN(size(sel));
        for i = 1:numel(sel)
            ZYX_sel(i) = ZYX(i,sel(i));
        end
        scatter(ax,per,log10(imag(ZYX_sel)),60,sel,'filled','MarkerEdgeColor','k','LineWidth',1.5)
    end 
    xlabel(ax,'period (s)')  
    ylabel(ax,'\Im log10 Zyx (mV/Km/nT)')
    title(ax,['\Im Zyx ',name(app.a)])

    set(ax,'Fontname','Verdana','Fontsize',12,...
           'box','on','Xscale','log','Yscale','linear',...
           'xlim',10.^[info.permin info.permax],...
           'xtick',10.^(ceil(info.permin):floor(info.permax)),...
           'MinorGridLineStyle','none')

    cb = colorbar(ax);
    cb.FontSize = 12;
    cb.Layout.Tile = 'east';
    cb.Label.String = 'Best time windows selected (%)';    
    cb.Label.FontName = 'Verdana';
    cb.Label.FontSize = 12;
    caxis(ax,[1,ranking])

end