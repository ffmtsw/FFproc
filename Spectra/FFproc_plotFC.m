function FFproc_plotFC(app,FC,tf,st,tw)
    
    global font 
    
    stat_vec = 1+(st-1)*5:5+(st-1)*5;
    if ~app.all.State
        FC_Original = cat(1,FC{tf}(1:end).FC_Original);
        FC_Mahalanobis = cat(1,FC{tf}(1:end).FC_Mahalanobis);
        FC_Fitted = cat(1,FC{tf}(1:end).FC_Fitted);
        RMS1 = cat(1,FC{tf}(1:end).RMS1);
        RMS2 = cat(1,FC{tf}(1:end).RMS2);
    else
        FC_Original = cat(1,FC{tf}(tw).FC_Original); 
        FC_Mahalanobis = cat(1,FC{tf}(tw).FC_Mahalanobis);
        FC_Fitted = cat(1,FC{tf}(tw).FC_Fitted);
        RMS1 = cat(1,FC{tf}(tw).RMS1);
        RMS2 = cat(1,FC{tf}(tw).RMS2);
    end
    
    FC_Original = FC_Original(:,stat_vec);
    FC_Mahalanobis = FC_Mahalanobis(:,stat_vec);
    FC_Fitted = FC_Fitted(:,stat_vec);
    RMS1 = RMS1(:,stat_vec);
    RMS2 = RMS2(:,stat_vec);
    
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
    
    %% EX
    ax = subplot(4,3,[1,4],'Parent',app.FCPanel);  
    cla(ax','reset')
    hold(ax,'on')
    grid(ax,'on')    
    
    if app.OriginalButton.Value
        scatter(ax,real(FC_Original(:,1)),imag(FC_Original(:,1)),4,'k','filled')    
    end
    if app.MahalanobisButton.Value
        scatter(ax,real(FC_Mahalanobis(:,1)),imag(FC_Mahalanobis(:,1)),4,'b','filled')
    end    
    if app.FittedButton.Value
        scatter(ax,real(FC_Fitted(:,1)),imag(FC_Fitted(:,1)),4,'r','filled')
    end  

    xlabel(ax,'\bf\Re EX (mV/Km)')  
    ylabel(ax,'\bf\Im EX (mV/Km)')
    set(ax,'Fontname',font,'Fontsize',10,...
           'box','on',...
           'XLim',[-max(abs(ax.XLim)) max(abs(ax.XLim))],...
           'YLim',[-max(abs(ax.YLim)) max(abs(ax.YLim))],...
           'MinorGridLineStyle','none')       
      
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
        
   %% EY
   ax = subplot(4,3,[7,10],'Parent',app.FCPanel);  
    cla(ax,'reset')
    hold(ax,'on')
    grid(ax,'on')
    
    if app.OriginalButton.Value
        scatter(ax,real(FC_Original(:,2)),imag(FC_Original(:,2)),4,'k','filled')        
    end
    if app.MahalanobisButton.Value
        scatter(ax,real(FC_Mahalanobis(:,2)),imag(FC_Mahalanobis(:,2)),4,'b','filled')
    end    
    if app.FittedButton.Value
        scatter(ax,real(FC_Fitted(:,2)),imag(FC_Fitted(:,2)),4,'r','filled')
    end        
        
    xlabel(ax,'\bf\Re EY (mV/Km)')  
    ylabel(ax,'\bf\Im EY (mV/Km)')
    set(ax,'Fontname',font,'Fontsize',10,...
           'box','on',...
           'XLim',[-max(abs(ax.XLim)) max(abs(ax.XLim))],...
           'YLim',[-max(abs(ax.YLim)) max(abs(ax.YLim))],...
           'MinorGridLineStyle','none')
       
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
       
    %% EY
    ax = subplot(4,3,[2,5],'Parent',app.FCPanel);  
    cla(ax,'reset')
    hold(ax,'on')
    grid(ax,'on')
    
    if app.OriginalButton.Value
        scatter(ax,real(FC_Original(:,3)),imag(FC_Original(:,3)),4,'k','filled')        
    end
    if app.MahalanobisButton.Value
        scatter(ax,real(FC_Mahalanobis(:,3)),imag(FC_Mahalanobis(:,3)),4,'b','filled')
    end    
    if app.FittedButton.Value
        scatter(ax,real(FC_Fitted(:,3)),imag(FC_Fitted(:,3)),4,'r','filled')
    end        
        
    xlabel(ax,'\bf\Re BX (nT)')  
    ylabel(ax,'\bf\Im BX (nT)')
    set(ax,'Fontname',font,'Fontsize',10,...
           'box','on',...
           'XLim',[-max(abs(ax.XLim)) max(abs(ax.XLim))],...
           'YLim',[-max(abs(ax.YLim)) max(abs(ax.YLim))],...
           'MinorGridLineStyle','none')
       
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
       
    %% BY
    ax = subplot(4,3,[8,11],'Parent',app.FCPanel);  
    cla(ax,'reset')
    hold(ax,'on')
    grid(ax,'on')
    
    if app.OriginalButton.Value
        scatter(ax,real(FC_Original(:,4)),imag(FC_Original(:,4)),4,'k','filled')        
    end
    if app.MahalanobisButton.Value
        scatter(ax,real(FC_Mahalanobis(:,4)),imag(FC_Mahalanobis(:,4)),4,'b','filled')
    end    
    if app.FittedButton.Value
        scatter(ax,real(FC_Fitted(:,4)),imag(FC_Fitted(:,4)),4,'r','filled')
    end        
        
    xlabel(ax,'\bf\Re BY (nT)')  
    ylabel(ax,'\bf\Im BY (nT)')
    set(ax,'Fontname',font,'Fontsize',10,...
           'box','on',...
           'XLim',[-max(abs(ax.XLim)) max(abs(ax.XLim))],...
           'YLim',[-max(abs(ax.YLim)) max(abs(ax.YLim))],...
           'MinorGridLineStyle','none')
       
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
       
    %% BZ
    ax = subplot(4,3,[3,6],'Parent',app.FCPanel);  
    cla(ax,'reset')
    hold(ax,'on')
    grid(ax,'on')
    
    if app.OriginalButton.Value
        scatter(ax,real(FC_Original(:,5)),imag(FC_Original(:,5)),4,'k','filled')        
    end
    if app.MahalanobisButton.Value
        scatter(ax,real(FC_Mahalanobis(:,5)),imag(FC_Mahalanobis(:,5)),4,'b','filled')
    end    
    if app.FittedButton.Value
        scatter(ax,real(FC_Fitted(:,5)),imag(FC_Fitted(:,5)),4,'r','filled')
    end        
        
    xlabel(ax,'\bf\Re BZ (nT)')  
    ylabel(ax,'\bf\Im BZ (nT)')
    set(ax,'Fontname',font,'Fontsize',10,...
           'box','on',...
           'XLim',[-max(abs(ax.XLim)) max(abs(ax.XLim))],...
           'YLim',[-max(abs(ax.YLim)) max(abs(ax.YLim))],...
           'MinorGridLineStyle','none')
       
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
    
    %% ERRORS RMS1
    ax = subplot(4,3,9,'Parent',app.FCPanel);  
    cla(ax,'reset')
    hold(ax,'on')
    grid(ax,'on')
    
    data = [median(RMS1,1);mad(RMS1,1,1);max(RMS1,[],1)];
    bar(ax, repmat(1:size(RMS1,2),3,1)',data')
    
    xlabel(ax,'\bfChannel')  
    ylabel(ax,'\bfRMSE_1')
    legend(ax,'Median','MAD','Max','FontName',font,'FontSize',8,'NumColumns',3);
    
    for i = 1:size(RMS1,2)/5
        chan{i,1} = ['Ex_',num2str(i)];
        chan{i,2} = ['Ey_',num2str(i)];
        chan{i,3} = ['Bx_',num2str(i)];
        chan{i,4} = ['By_',num2str(i)];
        chan{i,5} = ['Bz_',num2str(i)];
    end
    
    set(ax,'Fontname',font,'Fontsize',10,...
           'box','on',...
           'XLim',[0.5,numel(chan)+0.5],...
           'XTickLabels',chan(:),...
           'MinorGridLineStyle','none')  
    
    %% ERRORS RMS2
    ax = subplot(4,3,12,'Parent',app.FCPanel);  
    cla(ax,'reset')
    hold(ax,'on')
    grid(ax,'on') 
    
    data = [median(RMS2,1);mad(RMS2,1,1);max(RMS2,[],1)];
    bar(ax, repmat(1:size(RMS2,2),3,1)',data')
    
    xlabel(ax,'\bfChannel')  
    ylabel(ax,'\bfRMSE_2')
    
    set(ax,'Fontname',font,'Fontsize',10,...
           'box','on',...
           'XLim',[0.5,numel(chan)+0.5],...
           'XTickLabels',chan(:),...
           'MinorGridLineStyle','none')
    
end