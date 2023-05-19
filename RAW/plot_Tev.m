function plot_Tev(app,info,RAW,ranking,targ,id,cf)
    
    warning off 
    per = RAW(1).per;
    name = RAW(1).site;

    % Retrieving all estimations
    pertot = cat(1,RAW(1:ranking).per);
    TXZ = cat(1,RAW(1:ranking).txz);
    TYZ = cat(1,RAW(1:ranking).tyz);
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
    
    scatter(ax,pertot,real(TXZ(:,app.a)),50,ind,'filled') 
    if exist('id','var') && ~isempty(id)
        if ~isempty(cf(app.a).txz_Re)
            plot(ax,perdelta,cf(app.a).txz_Re(log10(perdelta)),'k:','LineWidth',1.5)
            for i = 1:numel(per)
                try
                    scatter(ax,per(i),real(RAW(id(app.a).txz_Re(i)).txz(i,app.a)),60,id(app.a).txz_Re(i),'filled','MarkerEdgeColor','k','LineWidth',1.5)
                catch
                end
            end
        end
    end   
    xlabel(ax,'period (s)')  
    ylabel(ax,'\Re Txz')
    title(ax,['\Re Txz ',name(app.a)])

    set(ax,'Fontname','Verdana','Fontsize',12,...
           'box','on','Xscale','log','Yscale','linear',...
           'xlim',10.^[info.permin info.permax],...
           'ylim',[-info.txz info.txz],...
           'xtick',10.^(ceil(info.permin):floor(info.permax)),...            
           'MinorGridLineStyle','none')
    
    % ------------------------------------------------------- Imaginary Zxy
    ax = nexttile(td);    
    cla(ax,'reset')
    hold(ax,'on')
    grid(ax,'on')
    colormap(ax,colors('viridis',ranking-1)) 
    
    scatter(ax,pertot,imag(TXZ(:,app.a)),50,ind,'filled') 
    if exist('id','var') && ~isempty(id)
        if ~isempty(cf(app.a).txz_Im)
            plot(ax,perdelta,cf(app.a).txz_Im(log10(perdelta)),'k:','LineWidth',1.5)
            for i = 1:numel(per)
                try
                    scatter(ax,per(i),imag(RAW(id(app.a).txz_Im(i)).txz(i,app.a)),60,id(app.a).txz_Im(i),'filled','MarkerEdgeColor','k','LineWidth',1.5)
                catch
                end
            end
        end
    end  
    xlabel(ax,'period (s)')  
    ylabel(ax,'\Im Txz')
    title(ax,['\Im Txz ',name(app.a)])

    set(ax,'Fontname','Verdana','Fontsize',12,...
           'box','on','Xscale','log','Yscale','linear',...
           'xlim',10.^[info.permin info.permax],...
           'ylim',[-info.txz info.txz],...
           'xtick',10.^(ceil(info.permin):floor(info.permax)),...            
           'MinorGridLineStyle','none')
    
    % ------------------------------------------------------------ Real Zxy
    ax = nexttile(td);     
    cla(ax,'reset')
    hold(ax,'on')
    grid(ax,'on')
    colormap(ax,colors('viridis',ranking-1))  
    
    scatter(ax,pertot,real(TYZ(:,app.a)),50,ind,'filled') 
    if exist('id','var') && ~isempty(id)
        if ~isempty(cf(app.a).tyz_Re)
            plot(ax,perdelta,cf(app.a).tyz_Re(log10(perdelta)),'k:','LineWidth',1.5)
            for i = 1:numel(per)
                try
                    scatter(ax,per(i),real(RAW(id(app.a).tyz_Re(i)).tyz(i,app.a)),60,id(app.a).tyz_Re(i),'filled','MarkerEdgeColor','k','LineWidth',1.5)
                catch
                end
            end
        end
    end 
    xlabel(ax,'period (s)')  
    ylabel(ax,'\Re Tyz')
    title(ax,['\Re Tyz ',name(app.a)])

    set(ax,'Fontname','Verdana','Fontsize',12,...
           'box','on','Xscale','log','Yscale','linear',...
           'xlim',10.^[info.permin info.permax],...
           'ylim',[-info.txz info.txz],...
           'xtick',10.^(ceil(info.permin):floor(info.permax)),...            
           'MinorGridLineStyle','none')
    
    % ------------------------------------------------------- Imaginary Zxy
    ax = nexttile(td);    
    cla(ax,'reset')
    hold(ax,'on')
    grid(ax,'on')    
    colormap(ax,colors('viridis',ranking-1)) 
    
    scatter(ax,pertot,imag(TYZ(:,app.a)),50,ind,'filled') 
    if exist('id','var') && ~isempty(id)
        if ~isempty(cf(app.a).tyz_Im)
            plot(ax,perdelta,cf(app.a).tyz_Im(log10(perdelta)),'k:','LineWidth',1.5)
            for i = 1:numel(per)
                try
                    scatter(ax,per(i),imag(RAW(id(app.a).tyz_Im(i)).tyz(i,app.a)),60,id(app.a).tyz_Im(i),'filled','MarkerEdgeColor','k','LineWidth',1.5)
                catch
                end
            end
        end
    end  
    xlabel(ax,'period (s)')  
    ylabel(ax,'\Im Tyz')
    title(ax,['\Im Txz ',name(app.a)])

    set(ax,'Fontname','Verdana','Fontsize',12,...
           'box','on','Xscale','log','Yscale','linear',...
           'xlim',10.^[info.permin info.permax],...
           'ylim',[-info.txz info.txz],...
           'xtick',10.^(ceil(info.permin):floor(info.permax)),...            
           'MinorGridLineStyle','none')
    
    cb = colorbar(ax);
    cb.FontSize = 12;
    cb.Layout.Tile = 'east';
    cb.Label.String = 'Best time windows selected (%)';    
    cb.Label.FontName = 'Verdana';
    cb.Label.FontSize = 12;
    caxis(ax,[0,ranking])

end