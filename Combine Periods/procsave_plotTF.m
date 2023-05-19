function procsave_plotTF(app,mt,info)
    
    font = 'Verdana';
    warning off
    % Indices from table
    tf = cat(1,app.table.Data{:,3});
    
    % Concatenating periods, resitivity and phase
    per = cat(1,mt.per);                % Period
    rhoxy = cat(1,mt.rhoxy);            % rhoxy
    rhoyx = cat(1,mt.rhoyx);            % rhoyx
    phixy = cat(1,mt.phixy);            % phixy
    phiyx = cat(1,mt.phiyx);            % phiyx
    txz = cat(1,mt.txz);                % txz
    tyz = cat(1,mt.tyz);                % tyz
    phimin = cat(1,mt.phimin);          % ptmin
    phimax = cat(1,mt.phimax);          % ptmax
    theta = cat(1,mt.theta);            % theta

    % Concatenating errors
    rhoxy_Err = cat(1,mt.rhoxy_Err);    % rhoxy_Err
    rhoyx_Err = cat(1,mt.rhoyx_Err);    % rhoyx_Err
    phixy_Err = cat(1,mt.phixy_Err);    % phixy_Err
    phiyx_Err = cat(1,mt.phiyx_Err);    % phiyx_Err
    txz_Err = cat(1,mt.txz_Err);        % txz_Err
    tyz_Err = cat(1,mt.tyz_Err);        % tyz_Err


    s1 = 18;
    s2 = 7;
    td = tiledlayout(2,50,'Parent',app.Panel,'TileSpacing','compact','Padding','compact');
      
    % -------------------------------------------------------- Rho diagonal
    ax1 = nexttile(td,[1,s1]);   
    cla(ax1)
    hold(ax1,'on')
    grid(ax1,'on')
    errorbar(ax1,per(tf),rhoxy(tf),...
                 rhoxy_Err(tf),...
                 'o','color',[0.7 0.1 0.1],'markersize',6,'markerfacecolor',[1,0.25,0.25],...
                 'capsize',3)
    errorbar(ax1,per(tf),rhoyx(tf),...
                 rhoyx_Err(tf),...
                 'o','color',[0.25,0.25,0.8],'markersize',6,'markerfacecolor',[0.1 0.6 0.9],...
                 'capsize',3)
    legend(ax1,'xy','yx','location','northwest','NumColumns',2,'Fontsize',8)
    xlabel(ax1,'period (s)')  
    ylabel(ax1,'\rho_a (\Omegam)')
    title(ax1,'Resistivity')

    set(ax1,'Fontname',font,'Fontsize',10,...
            'box','on','Xscale','log','Yscale','log',...
            'xlim',10.^[info.permin info.permax],...
            'ylim',10.^[info.rhomin info.rhomax],...         
            'MinorGridLineStyle','none')
    
    % ----------------------------------------------------------------- TXZ
    ax2 = nexttile(td,[1,s1]);   
    cla(ax2)
    hold(ax2,'on')
    grid(ax2,'on')
    errorbar(ax2,per(tf),real(txz(tf)),...
                 txz_Err(tf),...
                 'o','color',[0.7 0.1 0.1],'markersize',6,'markerfacecolor',[1,0.25,0.25],...
                 'capsize',3)
    errorbar(ax2,per(tf),imag(txz(tf)),...
                 txz_Err(tf),...
                 'o','color',[0.25,0.25,0.8],'markersize',6,'markerfacecolor',[0.1 0.6 0.9],...
                 'capsize',3)
    legend(ax2,'real','imag','location','northwest','NumColumns',2,'Fontsize',8)
    xlabel(ax2,'period (s)')  
    ylabel(ax2,'amplitude')
    title(ax2,'TXZ')

    set(ax2,'Fontname',font,'Fontsize',10,...
            'box','on','Xscale','log','Yscale','linear',...
            'xlim',10.^[info.permin info.permax],...
            'ylim',[info.tipmin info.tipmax],...           
            'MinorGridLineStyle','none')
        
    % -------------------------------------------------------------- Tipper
    ax5 = nexttile(td,[2,s2]);
    cla(ax5)
    hold(ax5,'on')
    grid(ax5,'on')
    
    quiver(ax5,0,log10(per(tf)),real(tyz(tf)),-real(txz(tf)),...
             'color',[1,0.25,0.25],'linewidth',1,'maxheadsize',0.1,'autoscale','off')
    
    quiver(ax5,0,log10(per(tf)),imag(tyz(tf)),-imag(txz(tf)),...
             'color',[0.1 0.6 0.9],'linewidth',1,'maxheadsize',0.1,'autoscale','off')
       
    set(ax5,'Fontname',font,'Fontsize',10,...
           'box','on',...
           'ydir','reverse',...
           'dataaspectratio',[1 1 1],...
           'xlim',[-0.6 0.6],'xtick',-5:0.5:5,'xticklabel',string(abs(-5:0.5:5)),...
           'ylim',[info.permin info.permax])
       
    xlabel(ax5,'amplitude')  
    ylabel(ax5,'log10 period (s)')
    title(ax5,'Tipper')
    legend(ax5,'Real','Imag','location','south','fontsize',8);

    % -------------------------------------------------------- Phase Tensor
    ax6 = nexttile(td,[2,s2]);
    
    cla(ax6,'reset')
    hold(ax6,'on')
    grid(ax6,'on')   
    
    % Limits for color scaling
    limmax = info.ptmax;
    limmin = info.ptmin;
        
    for i = 1:numel(tf)
        if tf(i)
            if 3*abs(atan(phimin(i))) < abs(atan(phimax(i)))
                ellmin = tan(atan(phimax(i))/3);
            else
                ellmin = phimin(i);
            end
            e = ellipse(0.1*atan(phimax(i))/atan(phimax(i)),...
                        0.1*atan(ellmin)/atan(phimax(i)),...
                        0,log10(per(i)),...
                        rad2deg(-theta(i)),...
                        'k',ax6);
            set(e,'facecolor','flat','facevertexcdata',rad2deg(atan(phimax(i))))
            
            b = minbar(0.9*0.1*atan(ellmin)/atan(phimax(i)),0.3*0.1,...
                       0,log10(per(i)),...
                       rad2deg(-theta(i)),...
                       'k',ax6);
            set(b,'facecolor','flat','facevertexcdata',rad2deg(atan(phimin(i))))  
        end
    end    
    
    set(ax6,'Fontname',font,'Fontsize',10,...
           'box','on',...
           'ydir','reverse',...
           'dataaspectratio',[1 1 1],...
           'xlim',[-0.6 0.6],'xtick',-5:0.5:5,...
           'ylim',[info.permin info.permax])
            
    title(ax6,'PT')
    ylabel(ax6,'log10 period (s)')
    
    clim(ax6,[limmin limmax])  
    
    colormap(ax6,colors('pt'))
    l = colorbar(ax6);
    l.Location = 'eastoutside';
    l.Label.String = '\Phi (°)';

    % ------------------------------------------------------ Phase diagonal
    ax3 = nexttile(td,[1,s1]);
    cla(ax3)
    hold(ax3,'on')
    grid(ax3,'on')
    
    errorbar(ax3,per(tf),phixy(tf),...
                 phixy_Err(tf),...
                 'o','color',[0.7 0.1 0.1],'markersize',6,'markerfacecolor',[1,0.25,0.25],...
                 'capsize',3)
    errorbar(ax3,per(tf),phiyx(tf)+180,...
                phiyx_Err(tf),...
                'o','color',[0.25,0.25,0.8],'markersize',6,'markerfacecolor',[0.1 0.6 0.9],...
                'capsize',3)
    legend(ax3,'xy','yx','location','northwest','NumColumns',2,'Fontsize',8)
    xlabel(ax3,'period (s)')  
    ylabel(ax3,'\phi (°)')
    title(ax3,'Phase')

    set(ax3,'Fontname',font,'Fontsize',10,...
            'box','on','Xscale','log','Yscale','linear',...
            'xlim',10.^[info.permin info.permax],...
            'ylim',[info.phimin info.phimax],...         
            'MinorGridLineStyle','none') 
        
    % ----------------------------------------------------------------- TYZ
    ax4 = nexttile(td,[1,s1]); 
    cla(ax4)
    hold(ax4,'on')
    grid(ax4,'on')
    errorbar(ax4,per(tf),real(tyz(tf)),...
                 tyz_Err(tf),...
                 'o','color',[0.7 0.1 0.1],'markersize',6,'markerfacecolor',[1,0.25,0.25],...
                 'capsize',3)
    errorbar(ax4,per(tf),imag(tyz(tf)),...
                 tyz_Err(tf),...
                 'o','color',[0.25,0.25,0.8],'markersize',6,'markerfacecolor',[0.1 0.6 0.9],...
                 'capsize',3)
    legend(ax4,'real','imag','location','northwest','NumColumns',2,'Fontsize',8)
    xlabel(ax4,'period (s)')  
    ylabel(ax4,'amplitude')
    title(ax4,'TYZ')

    set(ax4,'Fontname',font,'Fontsize',10,...
            'box','on','Xscale','log','Yscale','linear',...
            'xlim',10.^[info.permin info.permax],...
            'ylim',[info.tipmin info.tipmax],...         
            'MinorGridLineStyle','none')
        
    linkaxes([ax1,ax2,ax3,ax4],'x')    

end