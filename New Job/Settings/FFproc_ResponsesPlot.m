function FFproc_ResponsesPlot(app,ffts)
    
    fmax = app.maxfreqbox.Value;
    fmin = app.minfreqbox.Value;
    tfreq = mean(find_tfreq(app.targfreqbox.Value,fmin,fmax),1);
    for i = 1:numel(ffts)
        sr = ffts(i).sr;
        switch ffts(i).device
            case 'Metronix'                
                calibration = coil_cal_ADU(ffts(i).coil,sr);
                for j = 1:numel(ffts(i).coil)
                    amp{i,j} = interp1(calibration{j}(:,1),calibration{j}(:,2),tfreq','pchip');
                    phi{i,j} = interp1(calibration{j}(:,1),calibration{j}(:,3),tfreq','pchip');
                end
            case 'Phoenix'
                if sr == 24000
                    band = '2';
                elseif sr == 2400
                    band = '3';
                elseif sr == 150
                    band = '4';
                else
                    band = '5';
                end
                [~,filename,~] = fileparts([ffts(i).path,ffts(i).file]);
                ctsfile = [ffts(i).path,filename,'.CTS',band];
                [aresp,presp,~] = read_CTS(tfreq,ctsfile);
                for j = 1:numel(ffts(i).coil)
                    amp{i,j} = aresp(:,j+2);
                    phi{i,j} = rad2deg(presp(:,j+2));
                end
        end
        
    end
    
    fig = uifigure('Name','FFproc - Instrument Response','Position',[350,100,1000,850]);    
    tabgroup = uitabgroup(fig,'Position',[0,0,1000,850]);
    for i = 1:numel(ffts)
        % Bx Channel
        tab = uitab(tabgroup,'Title',['Coil - ',num2str(ffts(i).coil{1,1})]);
        tab.AutoResizeChildren = 'off';
        sp1 = subplot(2,1,1,'Parent',tab);
        plot(sp1,1./tfreq,amp{i,1},'o','color',[0.25,0.25,0.8],'markersize',6,'markerfacecolor',[0.1 0.6 0.9])
        title(sp1,['Site: ',ffts(i).name,' Bx-Channel'],'Interpreter','none')
        xlabel(sp1,'period (s)')
        set(sp1,'Xscale','Log','Fontname','Verdana','FontSize',12,...
                'xlim',10.^[floor(log10(min(1./tfreq))) ceil(log10(max(1./tfreq)))])
        grid(sp1,'on')
        sp2 = subplot(2,1,2,'Parent',tab);
        plot(sp2,1./tfreq,phi{i,1},'o','color',[0.25,0.25,0.8],'markersize',6,'markerfacecolor',[0.1 0.6 0.9])
        xlabel(sp2,'period (s)')
        ylabel(sp2,'phase \phi (°)')
        set(sp2,'Xscale','Log','Fontname','Verdana','FontSize',12,...
                'xlim',10.^[floor(log10(min(1./tfreq))) ceil(log10(max(1./tfreq)))],...
                'ylim',[0 90],...
                'ytick',0:15:90)
        grid(sp2,'on')

        % By Channel
        tab = uitab(tabgroup,'Title',['Coil - ',num2str(ffts(i).coil{1,2})]);
        tab.AutoResizeChildren = 'off';
        sp1 = subplot(2,1,1,'Parent',tab);
        plot(sp1,1./tfreq,amp{i,2},'o','color',[0.25,0.25,0.8],'markersize',6,'markerfacecolor',[0.1 0.6 0.9])
        title(sp1,['Site: ',ffts(i).name,' By-Channel'],'Interpreter','none')
        xlabel(sp1,'period (s)')
        set(sp1,'Xscale','Log','Fontname','Verdana','FontSize',12,...
                'xlim',10.^[floor(log10(min(1./tfreq))) ceil(log10(max(1./tfreq)))],...
                'ylim',10.^[floor(log10(min(min(cell2mat(amp(i,:)))))) ceil(log10(max(max(cell2mat(amp(i,:))))))])                
        grid(sp1,'on')
        sp2 = subplot(2,1,2,'Parent',tab);
        plot(sp2,1./tfreq,phi{i,2},'o','color',[0.25,0.25,0.8],'markersize',6,'markerfacecolor',[0.1 0.6 0.9])
        xlabel(sp2,'period (s)')
        ylabel(sp2,'phase \phi (°)')
        set(sp2,'Xscale','Log','Fontname','Verdana','FontSize',12,...
                'xlim',10.^[floor(log10(min(1./tfreq))) ceil(log10(max(1./tfreq)))],...
                'ylim',[0 90],...
                'ytick',0:15:90)
        grid(sp2,'on')

        % Bz Channel
        tab = uitab(tabgroup,'Title',['Coil - ',num2str(ffts(i).coil{1,3})]);
        tab.AutoResizeChildren = 'off';
        sp1 = subplot(2,1,1,'Parent',tab);
        plot(sp1,1./tfreq,amp{i,3},'o','color',[0.25,0.25,0.8],'markersize',6,'markerfacecolor',[0.1 0.6 0.9])
        title(sp1,['Site: ',ffts(i).name,' Bz-Channel'],'Interpreter','none')
        xlabel(sp1,'period (s)')
        set(sp1,'Xscale','Log','Fontname','Verdana','FontSize',12,...
                'xlim',10.^[floor(log10(min(1./tfreq))) ceil(log10(max(1./tfreq)))],...
                'ylim',10.^[floor(log10(min(min(cell2mat(amp(i,:)))))) ceil(log10(max(max(cell2mat(amp(i,:))))))])
        grid(sp1,'on')
        sp2 = subplot(2,1,2,'Parent',tab);
        plot(sp2,1./tfreq,phi{i,3},'o','color',[0.25,0.25,0.8],'markersize',6,'markerfacecolor',[0.1 0.6 0.9])
        xlabel(sp2,'period (s)')
        ylabel(sp2,'phase \phi (°)')
        set(sp2,'Xscale','Log','Fontname','Verdana','FontSize',12,...
                'xlim',10.^[floor(log10(min(1./tfreq))) ceil(log10(max(1./tfreq)))],...
                'ylim',[0 90],...
                'ytick',0:15:90)
        grid(sp2,'on')
    end
    

end
