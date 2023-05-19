% function read_CTS returns complex response functions from *.cts file
% created by Phoenix syscal_v6.exe for given set of frequencies

% version 1.1 / 14feb2019 / aj
% version 2.0 / 27may2020 / cc    Modified to fit in new version of FFproc,
%                                 name has been changed from MT_readCTS

% fr:       target frequencies (frequencies should be within range of cts
%           frequencies, otherwise NaN for response function
% ctsfile:  name of file containing Phoenix response functions
% resplim:  defines minimum value of response function amplitudes by
%           max(amplitude)*resplim

function resp = read_CTS(fr,ctsfile)

    cts = importdata(ctsfile,',',1);

    freq = cts.data(:,1);
    Nfreq = length(freq);
    chan = complex(cts.data(:,3:2:end),cts.data(:,4:2:end));
    
    if freq(end)-freq(1) < 0
        freq = flipud(freq);
        chan = flipud(chan);
    end

    Nchan = size(chan,2);
    achan = abs(chan);
    pchan = angle(chan)*180/pi;
    for nchan = 1:Nchan
        for nfreq = 2:Nfreq
            if pchan(nfreq,nchan)-pchan(nfreq-1,nchan) > 0
                pchan(nfreq:end,nchan) = pchan(nfreq:end,nchan) - 360;
            end
        end
    end

    aresp = interp1(freq,achan,fr.');
    presp = interp1(freq,pchan,fr.');
    resp = complex(aresp.*cosd(presp),aresp.*sind(presp));
    
end