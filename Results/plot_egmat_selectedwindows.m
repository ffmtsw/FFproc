function plot_egmat_selectedwindows(info,EV,targ,n)    
    
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
    
    for i = 1:numel(EV)
        vect_norm = 0:1/max(ntseg{i}):1;
        mat = zeros(max(ntseg{i}),length(tfreq{i}));
        win_thres = round(EV(i).perc.*EV(i).ntseg./100);
        for j = 1:length(tfreq{i})
            vec_ev = EV(i).timeseg(j).timeseg;
            mat_flag = zeros(size(vec_ev));
            mat_flag(vec_ev(1:win_thres(j))) = 1;
            mat_flag = [mat_flag(:);2];
            vec = 0:1/ntseg{i}(j):1;
            vec_int = interp1(vec,mat_flag,vect_norm,'previous')';
            mat(:,j) = vec_int(1:end-1,1);
        end
        MAT_flag{i} = mat;   
    end
    
    if isempty(targ)
        targ = figure('Name',['Captured Snapshots - ',num2str(EV(end).perc),'% Best Ranked Eigenvalues'],'WindowStyle','docked');
    end

    N = numel(extractfield(EV,'tfreq'));
    if numel(EV) == 1
        td = tiledlayout(1,1,'Parent',targ,'TileSpacing','tight','Padding','compact');
    else
        td = tiledlayout(N,1,'Parent',targ,'TileSpacing','compact','Padding','compact');
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

        imagesc(ax,vect{i},tfreq_v,MAT_flag{i}')

        if median(tfreq{i}) < 1
            set(ax,'YTick',(1:length(tfreq{i})),'YTickLabel',round(1./fliplr(tfreq{i}),2),...
                   'Ydir','reverse')
            ylabel(ax,'period (s)') 
        else
            set(ax,'YTick',(1:length(tfreq{i})),'YTickLabel',round(tfreq{i},2))
            ylabel(ax,'frequency (Hz)') 
        end
        title(ax,['EV Criteria ',' | ',num2str(EV(i).sr),' Hz | ',datestr(EV(i).tstart),' - ',datestr(EV(i).tend)])        
        xlim(ax,[min(vect{i}) max(vect{i})])
        ylim(ax,[min(tfreq_v)-0.5 max(tfreq_v)+0.5])

        % Colormap
        colormap(ax,flipud(bone(2)))
        % Colorbar
        if i == numel(EV)
            xlabel(ax,'time (h)')
            clb = colorbar(ax);
            clb.Layout.Tile = 'east';
            clb.Label.String = ['Time Window Selection by Eigenvalue Criteria - Best',num2str(EV(end).perc),'% Ranked Eigenvalues'];
            clb.Label.FontName = font;
            clb.Label.FontSize = fz;
            clb.Ticks = [0,1];
            clb.TickLabels = {'rejected','selected'};
        end

        caxis(ax,[0 1])
        set(ax,'fontname',font,'fontsize',fz,'box','on','layer','top')
    end     
    
end