function FFproc_TFInteractive(app,TF,tf,tw)
    
    global font 
    
    st1 = str2double(app.XStationTF.Value);
    st2 = str2double(app.YStationTF.Value);    
    
    XComp = reshape(extractfield(TF{tf},app.XTF.Value),numel(app.XStationTF.Items),[])';
    YComp = reshape(extractfield(TF{tf},app.YTF.Value),numel(app.YStationTF.Items),[])';
       
    if ~app.all.State
        switch app.XComponentTF.Value
            case 'Real'
                TF_X = real(XComp(:,st1));
                labelx = '\Re';
            case 'Imaginary'
                TF_X = imag(XComp(:,st1));
                labelx = '\Im';
        end
        switch app.YComponentTF.Value
            case 'Real'
                TF_Y = real(YComp(:,st2));
                labely = '\Re';
            case 'Imaginary'
                TF_Y = imag(XComp(:,st2));
                labely = '\Im';
        end
    else
        switch app.XComponentTF.Value
            case 'Real'
                TF_X = real(YComp(tw,st1));
                labelx = '\Re';
            case 'Imaginary'
                TF_X = imag(YComp(tw,st1));
                labelx = '\Im';
        end
        switch app.YComponentTF.Value
            case 'Real'
                TF_Y = real(YComp(tw,st2));
                labely = '\Re';
            case 'Imaginary'
                TF_Y = imag(YComp(tw,st2));
                labely = '\Im';
        end
    end 
        
    ax = subplot(1,1,1,'Parent',app.TFInteractivePanel);  
    cla(ax','reset')
    hold(ax,'on')
    grid(ax,'on')    

    scatter(ax,TF_X,TF_Y,15,1:numel(TF_X),'filled')  
    
    caxis(ax,[1,numel(app.TimeWindowSelectionListBox.Items)])
    colormap(ax,flipud(colors('lapaz')))
    cb = colorbar(ax);
    cb.Label.String = 'Time Window ID';
    
    xlabel(ax,['\bf',labelx,' ',app.XTF.Value,'_',num2str(st1)])  
    ylabel(ax,['\bf',labely,' ',app.YTF.Value,'_',num2str(st2)]) 
    set(ax,'Fontname',font,'Fontsize',14,...
           'box','on',...
           'XLim',[-max(abs(ax.XLim)) max(abs(ax.XLim))],...
           'YLim',[-max(abs(ax.YLim)) max(abs(ax.YLim))],...
           'MinorGridLineStyle','none') 

end