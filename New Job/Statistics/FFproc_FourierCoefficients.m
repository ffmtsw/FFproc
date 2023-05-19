clear, close, clc


% Univariate
[file,path] = uigetfile({'*.*', 'All Files (*.*)'; ...
                         '*.ts2;*.ts3;*.ts4;*.ts5;','Phoenix (*.ts*)';
                         '*.ats','Metronix (*.ats)';
                         '*.ts10','LEMI (*.t10)';
                         '*.RAW','EMERALD (*.RAW)';
                         '*.ML1;*.ML2;','Geolore (*.ML*)';
                         '*.asc','Uppsala ASCII (*.asc)';
                         '*.sec;*.min;','Intermagnet Observatory (*.sec,*.min)';
                         '*.mat','MatTS (*.mat)'},...
                         'Select input Time Series');     
% Load a new time series
ts = ts_input(path,file,true,[]);

% Time Series gaps
ts = findgaps(ts);   

% Number of oscillations
osc = 100;

% Target Frequencies
fmax = 6;
fmin = 0.001;
ftarg = 6;
tfreq_lim = find_tfreq(ftarg,fmin,fmax);
tfreq_vect = mean(tfreq_lim,2);

% Preallocating variables (vector)
n_tseg = NaN(numel(tfreq_vect),1);
n_samples = NaN(numel(tfreq_vect),1);
t_seg = NaN(numel(tfreq_vect),1);

%%
% clear FC
% Calculating spectrograms
FC = struct;
for f = 1:numel(tfreq_vect)

    clear F 
    F = cell(numel(ts),max(cat(1,ts(1).nchan)));
    
    [D,n_tseg(f),n_samples(f),t_seg(f),~,~,~,~] = ts_segment(ts,ts.sr,tfreq_lim(f,:),tfreq_vect(f),f,osc,false);

    % Samples id of frequency vector used to estimate each target frequency
    freq_id = round(tfreq_lim(f,:)*n_samples(f)/ts(1).sr); 

    % Calculating spectral lines for each target frequency
    % Frequency vector calculated as sr*(1:N)/N
    f_vec = ts(1).sr*(1:n_samples(f))/n_samples(f);
    
    % Segment of the frequency vector to be used
    freq_vec = f_vec(freq_id(2):freq_id(1)); 

    % Data calibration
    ...

    % Calculation of the Fourier Transform for each time segment (w), 
    % for each channel (k) and for each site (j)
    for j = 1:numel(ts)
        for k = 1:ts(j).nchan
            % Extracting each channel with n_tseg time segments
            dmat = D{j,k};     
            % Fourier Transform of each data channel Dmat
            Dmat = fft(dmat);
            % Creating calibration matrix Cmat
            Cmat = 1;
            % Matrix with calibrated Fourier Coefficients
            Fmat = Dmat.*Cmat;
            % Taking only 1/2 of the fft
            Fmat = Fmat(2:(end/2) + 1,:);
            % Cell array with j sites and k channels. Each entry
            % has only Fourier Coefficients around the target
            % frequency.
            F{j,k} = Fmat(freq_id(2):freq_id(1),:);
            % Releasing memory
            clear dmat Fmat Dmat Cmat
        end
        FC(f).EX(j).data = F{j,1};
        FC(f).EY(j).data = F{j,2};
        FC(f).BX(j).data = F{j,3};
        FC(f).BY(j).data = F{j,4};
        FC(f).BZ(j).data = F{j,5};
    end
end

%%
% Coherency between A and B
% coh(A,B) = <AB*><BA*>/<AA*><BB*>
COH = struct;
for j = 1:numel(ts)
    for f = 1:numel(tfreq_vect)    
        COH(f).EXBX(j).data = coh(FC(f).EX(j).data,FC(f).BX(j).data);
        COH(f).EYBX(j).data  = coh(FC(f).EY(j).data,FC(f).BX(j).data);
        COH(f).EXBY(j).data  = coh(FC(f).EX(j).data,FC(f).BY(j).data);
        COH(f).EYBY(j).data  = coh(FC(f).EY(j).data,FC(f).BY(j).data);
    end
end
%%
clear EXBX EXBY EYBX EYBY
maxL = max(n_tseg);
for j = 1:numel(ts)
    for f = 1:numel(tfreq_vect)    
        if numel(COH(f).EXBX(j)) <= 1
            EXBX{j}(f,:) = COH(f).EXBX(j).data;
            EYBX{j}(f,:) = COH(f).EYBX(j).data;
            EXBY{j}(f,:) = COH(f).EXBY(j).data;
            EXBY{j}(f,:) = COH(f).EXBY(j).data;
        else
            EXBX{j}(f,:) = interp1(COH(f).EXBX(j).data,linspace(1,n_tseg(f),maxL),'previous');
            EYBX{j}(f,:) = interp1(COH(f).EYBX(j).data,linspace(1,n_tseg(f),maxL),'previous');
            EYBX{j}(f,:) = interp1(COH(f).EYBX(j).data,linspace(1,n_tseg(f),maxL),'previous');
            EYBY{j}(f,:) = interp1(COH(f).EYBY(j).data,linspace(1,n_tseg(f),maxL),'previous');
        end    
    end
end

%%
tmin = ts.t(1);
tmax = ts.t(end);
t = linspace(tmin,tmax,maxL);

td = tiledlayout(2,2,"TileSpacing","tight");

for f = 1:numel(tfreq_vect)
    if tfreq_vect(f) < 1
        tfreqlabel{f} = [num2str(round(1./tfreq_vect(f),2)),' (s)'];
    else
        tfreqlabel{f} = [num2str(round(tfreq_vect(f),2)),' (Hz)'];
    end
end

[T,F] = meshgrid(t,1:numel(tfreq_vect));

nexttile(td)
pc = pcolor(T,F,EXBX);
pc.LineStyle = 'none';
ax = gca;
ylabel(ax,'coh^2 EXBX')
set(ax,'YDIr','reverse','Layer','Top','YTick',1:numel(tfreq_vect),'YTickLabel',tfreqlabel,'FontName','Verdana','FontSize',11)
colormap(ax,flipud(colors('viridis',40)))
caxis(ax,[0,1])

nexttile(td)
pc = pcolor(T,F,EYBX);
pc.LineStyle = 'none';
ax = gca;
ylabel(ax,'coh^2 EXBY')
set(ax,'YDIr','reverse','Layer','Top','XTickLabel',string(datetime(datevec(ax.XTick))),'YTick',1:numel(tfreq_vect),'YTickLabel',tfreqlabel,'FontName','Verdana','FontSize',11)
colormap(ax,flipud(colors('viridis',40)))
caxis(ax,[0,1])

nexttile(td)
pc = pcolor(T,F,EYBX);
pc.LineStyle = 'none';
ax = gca;
ylabel(ax,'coh^2 EYBX')
set(ax,'YDIr','reverse','Layer','Top','XTickLabel',string(datetime(datevec(ax.XTick))),'YTick',1:numel(tfreq_vect),'YTickLabel',tfreqlabel,'FontName','Verdana','FontSize',11)
colormap(ax,flipud(colors('viridis',40)))
caxis(ax,[0,1])

nexttile(td)
pc = pcolor(T,F,EYBY);
pc.LineStyle = 'none';
ax = gca;
ylabel(ax,'coh^2 EYBY')
set(ax,'YDIr','reverse','Layer','Top','XTickLabel',string(datetime(datevec(ax.XTick))),'YTick',1:numel(tfreq_vect),'YTickLabel',tfreqlabel,'FontName','Verdana','FontSize',11)
colormap(ax,flipud(colors('viridis',40)))
caxis(ax,[0,1])

cb = colorbar;
cb.Layout.Tile = 'east';
title(td,file,'FontName','Verdana','FontSize',16,'Interpreter','none')

%%
figure()
td = tiledlayout('flow',"TileSpacing","tight");
for f = 1:numel(tfreq_vect)
    nexttile(td)
    ax = gca;
%     t = linspace(tmin,tmax,);
%     tmed = [linspace(ts.t(1),ts())]
    scatter(datetime(datevec(linspace(tmin,tmax,n_tseg(f)))),COH(f).EXBX,5,'filled')
    title(ax,tfreqlabel(f),'FontName','Verdana','FontSize',11)
    set(ax,'FontName','Verdana','box','on')
    grid(ax,'on')
    xlim(ax,[datetime([datevec(tmin);datevec(tmax)])])
end
title(td,' \bfcoh^2 EXBX','FontName','Verdana','FontSize',14)

%%
figure()
td = tiledlayout('flow',"TileSpacing","tight");
for f = 1:numel(tfreq_vect)
    t = linspace(tmin,tmax,maxL);
    nexttile(td)
    ax = gca;
    scatter(1:n_tseg(f),COH(f).EXBY,5,'filled')
    title(ax,tfreqlabel(f),'FontName','Verdana','FontSize',11)
    set(ax,'FontName','Verdana','box','on')
end
title(td,' \bfcoh^2 EXBY','FontName','Verdana','FontSize',14)

%
figure()
td = tiledlayout('flow',"TileSpacing","tight");
for f = 1:numel(tfreq_vect)
    t = linspace(tmin,tmax,maxL);
    nexttile(td)
    ax = gca;
    scatter(1:n_tseg(f),COH(f).EYBX,5,'filled')
    title(ax,tfreqlabel(f),'FontName','Verdana','FontSize',11)
    set(ax,'FontName','Verdana','box','on')
end
title(td,' \bfcoh^2 EYBX','FontName','Verdana','FontSize',14)

%
figure()
td = tiledlayout('flow',"TileSpacing","tight");
for f = 1:numel(tfreq_vect)
    t = linspace(tmin,tmax,maxL);
    nexttile(td)
    ax = gca;
    scatter(1:n_tseg(f),COH(f).EYBY,5,'filled')
    title(ax,tfreqlabel(f),'FontName','Verdana','FontSize',11)
    set(ax,'FontName','Verdana','box','on')
end
title(td,' \bfcoh^2 EYBY','FontName','Verdana','FontSize',14)

%%
figure()
td = tiledlayout('flow',"TileSpacing","tight");
for f = 1:numel(tfreq_vect)
    nexttile(td)
    ax = gca;
%     t = linspace(tmin,tmax,);
%     tmed = [linspace(ts.t(1),ts())]
    scatter(datetime(datevec(linspace(tmin,tmax,n_tseg(f)))),COH(f).EXBX,5,'filled')
    title(ax,tfreqlabel(f),'FontName','Verdana','FontSize',11)
    set(ax,'FontName','Verdana','box','on')
    grid(ax,'on')
    xlim(ax,[datetime([datevec(tmin);datevec(tmax)])])
end
title(td,' \bfcoh^2 EXBX','FontName','Verdana','FontSize',14)

%%
function C = coh(A,B)

    % coh2(A,B) = <AB*><BA*>/<AA*><BB*>
    ABconj = sum(A.*conj(B),1);
    BAconj = sum(B.*conj(A),1);
    AAconj = sum(A.*conj(A),1);
    BBconj = sum(B.*conj(B),1);
    C = abs((ABconj.*BAconj)./(AAconj.*BBconj));

end

% Bivariate coherence between inputs B1 and B2 and output A
function bvC = bivarcoh(A,B1,B2)
    
    % bvcoh = (<B1A*> + <B2A*>)/<AA*>
    B1Aconj = sum(B1.*conj(A),1);
    B2Aconj = sum(B2.*conj(A),1);
    AAconj = sum(A.*conj(A),1);
    bvC = (B1Aconj + B2Aconj)./AAconj;

end