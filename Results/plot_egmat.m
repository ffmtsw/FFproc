function plot_egmat(info,EV,targ,n)
    
    global font 
    
    [~,ind] = sort(cat(1,EV.sr),'descend');   
    EV = EV(ind);
    fz = 11;
    
    for i = 1:numel(EV)
        tseg{i} = EV(i).tseg;
        ntseg{i} = EV(i).ntseg;
        tfreq{i} = EV(i).tfreq; 
        vect{i} = linspace(0,max((ntseg{i}+1).*tseg{i})/3600/2,max(ntseg{i}));
    end
    
    vect_norm = 0:1/999:1;
    
    for i = 1:numel(EV)
        mat = zeros(1000,length(tfreq{i}));
        for j = 1:length(tfreq{i})
            vec_ev = EV(i).EVal(j).EVal(:,n);
            vec_ev = [vec_ev;1];
            vec = 0:1/ntseg{i}(j):1;
            vec_int = interp1(vec,vec_ev,vect_norm,'previous');
            mat(:,j) = vec_int;
        end
        MAT_EV{i} = fliplr(mat);     
    end
    
    best = cat(1,EV.best);
    M = real(floor(log10(max(max(best)))));
%     m = log10(1);
    m = real(ceil(log10(min(min(best)))));
    
    for i = 1:numel(EV)
    % Upper Limit
        if max(max(abs(MAT_EV{i})))>10
            cmin(i) = 1;
        else
            cmin(i) = min(min(abs(MAT_EV{i})));
        end
        % Lower Limit
        if max(max(abs(MAT_EV{i})))<100
            cmax(i) = 100;
        else
            cmax(i) = max(max(abs(MAT_EV{i})));
        end
    end

    N = numel(extractfield(EV,'tfreq'));
    if numel(EV) == 1
        td = tiledlayout(1,1,'Parent',targ,'TileSpacing','tight','Padding','compact');
    else
        td = tiledlayout(N,1,'Parent',targ,'TileSpacing','tight','Padding','compact');
    end

    for i = 1:numel(EV)
        
        if numel(EV) == 1
            ax = nexttile(td);
        else
            ax = nexttile(td,[numel(EV(i).tfreq),1]); 
        end
        cla(ax,'reset')
        hold(ax,'on')
        grid(ax,'on')        
        
        if median(tfreq{i}) < 1
        tfreq_v = fliplr(1:length(tfreq{i}));
        else
            tfreq_v = (1:length(tfreq{i}));
        end

        imagesc(ax,vect{i},tfreq_v,abs(MAT_EV{i})')

        if median(tfreq{i}) < 1
            set(ax,'YTick',(1:length(tfreq{i})),'YTickLabel',round(1./fliplr(tfreq{i}),2),'Ydir','reverse')
            ylabel(ax,'period (s)') 
        else
            set(ax,'YTick',(1:length(tfreq{i})),'YTickLabel',round(tfreq{i},2),'Ydir','normal')
            ylabel(ax,'frequency (Hz)') 
        end
        title(ax,['Eigenvalue ',num2str(n),' | ',num2str(EV(i).sr),' Hz | ',datestr(EV(i).tstart),' - ',datestr(EV(i).tend)])
        xlim(ax,[min(vect{i}) max(vect{i})])
        ylim(ax,[min(tfreq_v)-0.5 max(tfreq_v)+0.5])

        % Colormap
        colormap(ax,colors('pastel',31))

        % Colorbar
        if i == numel(EV)
            xlabel(ax,'time (h)')
            clb = colorbar(ax);
            clb.Layout.Tile = 'east';
            clb.Label.String = 'Signal / Noise';
            clb.Label.FontName = font;
            clb.Label.FontSize = fz;
        end

        caxis(ax,[10^m 10^M])
        set(ax,'ColorScale','log','fontname',font,'fontsize',fz,...
                                             'box','on','layer','top')
    end     
    
end