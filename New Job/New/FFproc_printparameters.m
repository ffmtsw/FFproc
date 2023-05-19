function txt = FFproc_printparameters(job)

    % Extracting information from final structure
    target_freq = [job.minf,job.maxf];
    
    % Number of time series per job
    nts = numel(job.ffts);
    
    % Magnetic reference channels
    magref = zeros(numel(job.ffts),2);
    for j = 1:numel(job.ffts)
        magref(j,1) = job.ffts(j).rsbx;
        magref(j,2) = job.ffts(j).rsby;
    end
        
    
    %% Printing information
    txt{1,1} = '************************ TIME SERIES ************************';    
    txt{end+1,1} = '';
    txt{end+1,1} = ['Number of time series:   ',num2str(nts)];
    txt{end+1,1} = '';
        for i = 1:nts
            txt{end+1,1} = ['------------------------ Time series ',num2str(i),' ----------------------'];
            if iscell(job.ffts(i).file)
                for j = 1:numel(job.ffts(i).file)
                    txt{end+1,1} = ['Name:                ',job.ffts(i).file{j}];
                end
            else
                txt{end+1,1} = ['Name:                ',job.ffts(i).file];                
            end
            txt{end+1,1} = ['Site:                ',job.ffts(i).name];
            txt{end+1,1} = ['Directory:           ',job.ffts(i).path];
            txt{end+1,1} = '';
            txt{end+1,1} = ['Sampling rate (Hz):  ',num2str(job.sr)];
            txt{end+1,1} = ['Device:              ',job.ffts(i).device];
            txt{end+1,1} = ['Start date/time UTC: ',datestr(job.ffts(i).tstart)];
            txt{end+1,1} = ['End date/time UTC:   ',datestr(job.ffts(i).tend)];
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
            txt{end+1,1} = ['Channels:            ',ex,ey,bx,by,bz];
            if all(magref(i,:) == i)
                rr = 'false';
            else
                rr = 'true';
            end
            txt{end+1,1} = ['Remote Reference:    ',rr];
            txt{end+1,1} = ['Ex dipole (m):       ',num2str(job.ffts(i).diplen(1))];
            txt{end+1,1} = ['Ey dipole (m):       ',num2str(job.ffts(i).diplen(2))];
            txt{end+1,1} = ['E-azimuth (°):       ',num2str(job.ffts(i).rot)];
            if iscell(job.ffts(i).coil)
                try 
                    coil = cell2mat(job.ffts(i).coil);
                    txt{end+1,1} = ['Bx SN:               ',num2str(coil(1))];
                    txt{end+1,1} = ['By SN:               ',num2str(coil(2))];
                    txt{end+1,1} = ['Bz SN:               ',num2str(coil(3))]; 
                catch
                    txt{end+1,1} = ['Bx SN:               ',job.ffts(i).coil{1}];
                    txt{end+1,1} = ['By SN:               ',job.ffts(i).coil{2}];
                    txt{end+1,1} = ['Bz SN:               ',job.ffts(i).coil{3}];    
                end
            else
                txt{end+1,1} = ['Bx SN:               ',num2str(job.ffts(i).coil(1))];
                txt{end+1,1} = ['By SN:               ',num2str(job.ffts(i).coil(2))];
                txt{end+1,1} = ['Bz SN:               ',num2str(job.ffts(i).coil(3))];                
            end
            
            if job.ffts(i).RAP
                rap = 'true';
            else
                rap = 'false';
            end
            txt{end+1,1} = ['RAP:                 ',rap];   
            txt{end+1,1} = ['Latitude (°):        ',num2str(job.ffts(i).GPS(1))];
            txt{end+1,1} = ['Longitude (°):       ',num2str(job.ffts(i).GPS(2))];
            txt{end+1,1} = ['Elevation (m):       ',num2str(job.ffts(i).GPS(3))];
            txt{end+1,1} = '';
        end
        txt{end+1,1} = '';
        txt{end+1,1} = '************************* PARAMETERS ************************';
        txt{end+1,1} = '';
        txt{end+1,1} = '---------------------- Target Frequencies -------------------';
        txt{end+1,1} = ['Maximum Frequency (Hz):     ',num2str(job.maxf)];
        txt{end+1,1} = ['Minimum Frequency (Hz):     ',num2str(job.minf)];
        txt{end+1,1} = ['Frequencies per decade:     ',num2str(job.ftarg)];
        txt{end+1,1} = '';
        txt{end+1,1} = '-------------------------- Filtering ------------------------';
        txt{end+1,1} = ['Single Frequencies (Hz):    ',num2str(job.singf)];
        txt{end+1,1} = ['Notch Filter Width (Hz):    ',num2str(job.phsing)];
        txt{end+1,1} = ['Multiples (Hz):             ',num2str(job.mult)];
        txt{end+1,1} = ['Notch Filter Width (Hz):    ',num2str(job.phmult)];
        txt{end+1,1} = '';
        txt{end+1,1} = '------------------------ Oscillations -----------------------';
        txt{end+1,1} = ['Oscillations:               ',num2str(job.osc)];
        txt{end+1,1} = '';
        txt{end+1,1} = '------------------------- Eigenvalues -----------------------';
        txt{end+1,1} = ['EV Index Criteria:          ',num2str(job.eind)];
        txt{end+1,1} = ['EV Criteria (%):            ',num2str(job.ranking)];
        txt{end+1,1} = '';
        txt{end+1,1} = '---------------------- Robust Parameters --------------------';
        txt{end+1,1} = ['Regresion:                  ',job.reg];
        txt{end+1,1} = ['Noise Model Estimator:      ',job.model];
        txt{end+1,1} = '';
        txt{end+1,1} = '*************************** OUTPUT **************************';
        txt{end+1,1} = '';
        txt{end+1,1} = 'OPT:                        true';
        if job.rnk
            rnk = 'true';
        else
            rnk = 'false';
        end
        txt{end+1,1} = ['RNK:                        ',rnk];
        if job.raw
            raw = 'true';
        else
            raw = 'false';
        end
        txt{end+1,1} = ['RAW:                        ',raw];
        if job.fc
            fc = 'true';
        else
            fc = 'false';
        end
        txt{end+1,1} = ['FC:                         ',fc];

end