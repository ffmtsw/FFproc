function FFproc_FCInteractive(app,FC,tf,tw)
    
    global font 
    
    st1 = str2double(app.XStation.Value);
    st2 = str2double(app.YStation.Value);
    
    switch app.XChannel.Value
        case 'Ex'
            ch1 = 1+(st1-1)*5;
            unit1 = ' (mV/Km)';
        case 'Ey'
            ch1 = 2+(st1-1)*5;
            unit1 = ' (mV/Km)';
        case 'Bx'
            ch1 = 3+(st1-1)*5;
            unit1 = ' (nT)';
        case 'By'
            ch1 = 4+(st1-1)*5;
            unit1 = ' (nT)';
        case 'Bz'
            ch1 = 5+(st1-1)*5;
            unit1 = ' (nT)';
    end
    switch app.YChannel.Value
        case 'Ex'
            ch2 = 1+(st2-1)*5;
            unit2 = ' (mV/Km)';
        case 'Ey'
            ch2 = 2+(st2-1)*5;
            unit2 = ' (mV/Km)';
        case 'Bx'
            ch2 = 3+(st2-1)*5;
            unit2 = ' (nT)';
        case 'By'
            ch2 = 4+(st2-1)*5;
            unit2 = ' (nT)';
        case 'Bz'
            ch2 = 5+(st2-1)*5;
            unit2 = ' (nT)';
    end    
       
    if ~app.all.State
        switch app.XComponent.Value
            case 'Real'
                FC_Original_1 = real(cat(1,FC{tf}(1:end).FC_Original));
                FC_Mahalanobis_1 = real(cat(1,FC{tf}(1:end).FC_Mahalanobis));
                FC_Fitted_1 = real(cat(1,FC{tf}(1:end).FC_Fitted));
                labelx = '\Re';
            case 'Imaginary'
                FC_Original_1 = imag(cat(1,FC{tf}(1:end).FC_Original));
                FC_Mahalanobis_1 = imag(cat(1,FC{tf}(1:end).FC_Mahalanobis));
                FC_Fitted_1 = imag(cat(1,FC{tf}(1:end).FC_Fitted)); 
                labelx = '\Im';
        end
        switch app.YComponent.Value
            case 'Real'
                FC_Original_2 = real(cat(1,FC{tf}(1:end).FC_Original));
                FC_Mahalanobis_2 = real(cat(1,FC{tf}(1:end).FC_Mahalanobis));
                FC_Fitted_2 = real(cat(1,FC{tf}(1:end).FC_Fitted));
                labely = '\Re';
            case 'Imaginary'
                FC_Original_2 = imag(cat(1,FC{tf}(1:end).FC_Original));
                FC_Mahalanobis_2 = imag(cat(1,FC{tf}(1:end).FC_Mahalanobis));
                FC_Fitted_2 = imag(cat(1,FC{tf}(1:end).FC_Fitted)); 
                labely = '\Im';
        end
    else
        switch app.XComponent.Value
            case 'Real'
                FC_Original_1 = cat(1,FC{tf}(tw).FC_Original);
                FC_Mahalanobis_1 = real(cat(1,FC{tf}(tw).FC_Mahalanobis));
                FC_Fitted_1 = real(cat(1,FC{tf}(tw).FC_Fitted));
                labelx = '\Re';
            case 'Imaginary'
                FC_Original_1 = imag(cat(1,FC{tf}(tw).FC_Original));
                FC_Mahalanobis_1 = imag(cat(1,FC{tf}(tw).FC_Mahalanobis));
                FC_Fitted_1 = imag(cat(1,FC{tf}(tw).FC_Fitted));
                labelx = '\Im';
        end
        switch app.YComponent.Value
            case 'Real'
                FC_Original_2 = real(cat(1,FC{tf}(tw).FC_Original));
                FC_Mahalanobis_2 = real(cat(1,FC{tf}(tw).FC_Mahalanobis));
                FC_Fitted_2 = real(cat(1,FC{tf}(tw).FC_Fitted));
                labely = '\Re';
            case 'Imaginary'
                FC_Original_2 = imag(cat(1,FC{tf}(tw).FC_Original));
                FC_Mahalanobis_2 = imag(cat(1,FC{tf}(tw).FC_Mahalanobis));
                FC_Fitted_2 = imag(cat(1,FC{tf}(tw).FC_Fitted)); 
                labely = '\Im';
        end
    end 
    
    FC_Original_1 = FC_Original_1(:,ch1);
    FC_Mahalanobis_1 = FC_Mahalanobis_1(:,ch1);
    FC_Fitted_1 = FC_Fitted_1(:,ch1);
    FC_Original_2 = FC_Original_2(:,ch2);
    FC_Mahalanobis_2 = FC_Mahalanobis_2(:,ch2);
    FC_Fitted_2 = FC_Fitted_2(:,ch2);
    
    ax = subplot(1,1,1,'Parent',app.FCInteractivePanel);  
    cla(ax','reset')
    hold(ax,'on')
    grid(ax,'on')    
    if app.OriginalButton.Value
        scatter(ax,FC_Original_1,FC_Original_2,6,'k','filled')    
    end
    if app.MahalanobisButton.Value
        scatter(ax,FC_Mahalanobis_1,FC_Mahalanobis_2,5,'b','filled')
    end    
    if app.FittedButton.Value
        scatter(ax,FC_Fitted_1,FC_Fitted_2,4,'r','filled')
    end  

    xlabel(ax,['\bf',labelx,' ',app.XChannel.Value,'_',num2str(st1),unit1])  
    ylabel(ax,['\bf',labely,' ',app.YChannel.Value,'_',num2str(st2),unit2]) 
    set(ax,'Fontname',font,'Fontsize',14,...
           'box','on',...
           'XLim',[-max(abs(ax.XLim)) max(abs(ax.XLim))],...
           'YLim',[-max(abs(ax.YLim)) max(abs(ax.YLim))],...
           'MinorGridLineStyle','none') 
    
    if app.OriginalButton.Value
        original = true;
    else
        original = false;
    end
    if app.MahalanobisButton.Value
        maha = true;
    else
        maha = false;
    end    
    if app.FittedButton.Value
        fitted = true;
    else
        fitted = false;
    end 
    if original && maha && fitted
        leg = legend(ax,'Original Signal','Mahalanobis OD','Estimated Signal','FontName',font,'FontSize',8);
    elseif original && maha
        leg = legend(ax,'Original Signal','Mahalanobis OD','FontName',font,'FontSize',8);
    elseif maha && fitted
        leg = legend(ax,'Mahalanobis OD','Estimated Signal','FontName',font,'FontSize',8);
    elseif original && fitted
        leg = legend(ax,'Original Signal','Estimated Signal','FontName',font,'FontSize',8);
    elseif original
        leg = legend(ax,'Original Signal','FontName',font,'FontSize',8);
    elseif maha
        leg = legend(ax,'Mahalanobis OD','FontName',font,'FontSize',8);
    elseif fitted
        leg = legend(ax,'Estimated Signal','FontName',font,'FontSize',8);
    else
    end
end