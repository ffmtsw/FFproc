function FFproc_PRM_save(final)

    % Extracting information from final structure
    job = final.job;
    % Target frequencies
    target_freq = final.target_freq;
    % Number of oscillations
    osc = job.osc;

    % Processing approach
    if job.bivar && ~job.ev         % Bivariate
        approach = 'BVR';
    elseif ~job.bivar && job.ev     % Eigenvalue Decomposition
        approach = 'EVD';
        approach = [approach,'-',job.eind];
    end

    ranking = num2str(job.ranking);
    if numel(ranking) == 1
        ranking = ['00',ranking];
    elseif numel(ranking) == 2
        ranking = ['0',ranking];
    end

    % Transfer functions ranking
    zrank = ['ZRNK-',job.zrank];
    trank = ['TRNK-',job.trank];

    % Powers
    evpow = ['EIPOW-',num2str(job.evpow)];
    pcpow = ['PCPOW-',num2str(job.rpow)];
    
    % Number of time series per job
    nts = numel(job.ffts);
    
    % Magnetic reference channels
    magref = zeros(numel(job.ffts),2);
    for j = 1:numel(job.ffts)
        magref(j,1) = job.ffts(j).rsbx;
        magref(j,2) = job.ffts(j).rsby;
    end
    
    %% Saving estimations in corresponding files
    load('ffprocsavepath.mat','path')
    cd(path)
    
    % Creating folder name for each job an creating folder (if it does not
    % exists)
    file = job.name;
    foldername = [];
    if numel(file) == 1
        foldername = file{1};
    else
        for n = 1:numel(job.name)
            foldername = [foldername,[file{n}]]; %#ok<*AGROW>
            if n < numel(file)
                foldername = [foldername,'-'];
            end
        end
    end
    if isfolder(foldername)
    else
        mkdir(foldername)
    end
    cd(foldername)
    
    % Retrieving target frequency range label: seconds or Hz
    if max(target_freq) >= 1
        f1 = round(target_freq(1));
        unit1 = 'Hz';
    else
        f1 = round(1./target_freq(1));
        unit1 = 's';
    end
    if min(target_freq) >= 1
        f2 = round(target_freq(end));
        unit2 = 'Hz';
    else
        f2 = round(1./target_freq(end));
        unit2 = 's';
    end
    
    % Remote Reference indicator
    if  numel(magref(:,1)) ~= numel(unique(magref(:,1))) || ...
        numel(magref(:,2)) ~= numel(unique(magref(:,2)))
        remote = '_RR';
    else
        remote = '';
    end
    
    switch job.model
        case 'Standard'
            model = 'STD_';
        case 'Advanced'
            model = 'ANM_';
    end

    switch job.reg
        case 'Robust'
            reg = 'ROB_';
        case 'LSQ'
            reg = 'LSQ_';
    end
    
    % Assembling parameters file name
    savename = [model,reg,num2str(osc),'OSC_',...
                num2str(f1),unit1,'-',num2str(f2),unit2,...
                remote,'_'...
                approach,'_',zrank,'_',trank,'_',evpow,'_',pcpow,...
                '_PRM',ranking,'%.txt'];
    
    %% Printing information
    fid = fopen(savename,'wt');
    
    try
        fprintf(fid,'************************ TIME SERIES ************************\n\n');    
        fprintf(fid,['Number of time series:\t',num2str(nts),'\n']);

        for i = 1:nts
            fprintf(fid,['\n------------------------ Time series ',num2str(i),' ----------------------\n']);
            if iscell(job.ffts(i).file)
                for j = 1:numel(job.ffts(i).file)
                    fprintf(fid,['Name:                ',job.ffts(i).file{j}]);
                end
                fprintf(fid,'\n');
            else
                fprintf(fid,['Name:                ',job.ffts(i).file,'\n']);                
            end
            fprintf(fid,['Site:                ',job.ffts(i).name,'\n']);
            fprintf(fid,'%s',['Directory:           ',job.ffts(i).path]);
            fprintf(fid,'\n');
            fprintf(fid,['Sampling rate (Hz):  ',num2str(job.sr),'\n']);
            fprintf(fid,['Device:              ',job.ffts(i).device,'\n']);
            fprintf(fid,['Start date/time UTC: ',datestr(job.ffts(i).tstart),'\n']);
            fprintf(fid,['End date/time UTC:   ',datestr(job.ffts(i).tend),'\n']);
            if job.ffts(i).actchan(1)
                ex = 'Ex ';
            else
                ex = '';
            end
            if job.ffts(i).actchan(2)
                ey = 'Ey ';
            else
                ey = '';
            end
            if job.ffts(i).actchan(3)
                bx = 'Bx ';
            else
                bx = '';
            end
            if job.ffts(i).actchan(4)
                by = 'By ';
            else
                by = '';
            end
            if job.ffts(i).actchan(5)
                bz = 'Bz ';
            else
                bz = '';
            end
            fprintf(fid,['Channels:            ',ex,ey,bx,by,bz,'\n']);
            if all(magref(i,:) == i)
                rr = 'false';
            else
                rr = 'true';
            end
            fprintf(fid,['Remote Reference:    ',rr,'\n']);
            fprintf(fid,['Ex dipole (m):       ',num2str(job.ffts(i).diplen(1)),'\n']);
            fprintf(fid,['Ey dipole (m):       ',num2str(job.ffts(i).diplen(2)),'\n']);
            fprintf(fid,['E-azimuth (°):       ',num2str(job.ffts(i).rot),'\n']);
            fprintf(fid,['Bx SN:               ',job.ffts(1).coil{1},'\n']);
            fprintf(fid,['By SN:               ',job.ffts(1).coil{2},'\n']);
            fprintf(fid,['Bz SN:               ',job.ffts(1).coil{3},'\n']);
            if job.ffts(i).RAP
                rap = 'true';
            else
                rap = 'false';
            end
            fprintf(fid,['RAP:                 ',rap,'\n']);   
            fprintf(fid,['Latitude (°):        ',num2str(job.ffts(i).GPS(1)),'\n']);
            fprintf(fid,['Longitude (°):       ',num2str(job.ffts(i).GPS(2)),'\n']);
            fprintf(fid,['Elevation (m):       ',num2str(job.ffts(i).GPS(3)),'\n']);
        end
        fprintf(fid,'\n');
        fprintf(fid,'\n************************* PARAMETERS ************************\n');
        fprintf(fid,'\n---------------------- Target Frequencies -------------------\n');
        fprintf(fid,['Maximum Frequency (Hz):     ',num2str(job.maxf),'\n']);
        fprintf(fid,['Minimum Frequency (Hz):     ',num2str(job.minf),'\n']);
        fprintf(fid,['Frequencies per decade:     ',num2str(job.ftarg),'\n']);
        fprintf(fid,'\n');
        fprintf(fid,'\n-------------------------- Filtering ------------------------\n');
        fprintf(fid,['Single Frequencies (Hz):    ',num2str(job.singf),'\n']);
        fprintf(fid,['Notch Filter Width (Hz):    ',num2str(job.phsing),'\n']);
        fprintf(fid,['Multiples (Hz):             ',num2str(job.mult),'\n']);
        fprintf(fid,['Notch Filter Width (Hz):    ',num2str(job.phmult),'\n']);
        fprintf(fid,'\n');
        fprintf(fid,'\n------------------------ Oscillations -----------------------\n');
        fprintf(fid,['Oscillations:               ',num2str(job.osc),'\n']);
        fprintf(fid,'\n');
        fprintf(fid,'\n------------------------- Eigenvalues -----------------------\n');
        fprintf(fid,['EV Index Criteria:          ',num2str(job.eind),'\n']);
        fprintf(fid,['EV Criteria (%):            ',num2str(job.ranking),'\n']);
        fprintf(fid,'\n');
        fprintf(fid,'\n---------------------- Robust Parameters --------------------\n');
        fprintf(fid,['Regresion:                  ',job.reg,'\n']);
        fprintf(fid,['Noise Model Estimator:      ',job.model,'\n']);
        fprintf(fid,'\n');
        fprintf(fid,'\n*************************** OUTPUT **************************\n');
        fprintf(fid,['OPT:                        true','\n']);
        if job.raw
            raw = 'true';
        else
            raw = 'false';
        end
        fprintf(fid,['RAW:                        ',raw,'\n']);
        if job.fc
            fc = 'true';
        else
            fc = 'false';
        end
        fprintf(fid,['FC:                         ',fc,'\n']);
        fprintf(fid,'%s',['Directory:                  ',pwd]);   
        fclose(fid);
    catch
        fclose(fid);
    end

end