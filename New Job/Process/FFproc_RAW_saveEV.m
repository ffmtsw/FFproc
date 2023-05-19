function EV = FFproc_RAW_saveEV(job,target_freq,ntseg,tseg,Results,EIGVAL)  

    for i = 1:numel(Results)
        EV.best(i).best = Results(i).EVal;
    end
    
    EV.perc = job.eigperc;
    EV.sr = job.sr;
    
    for i = 1:numel(target_freq)
        EV.EVal(i).EVal = EIGVAL(i).EVal;
        EV.EV(i).EV = EIGVAL(i).EV;
        EV.timeseg(i).timeseg = EIGVAL(i).timeseg;
    end
    
    EV.tfreq = target_freq;
    EV.ntseg = ntseg;
    EV.tseg = tseg;
    EV.tstart = job.ffts(1).tstart;
    EV.tend = job.ffts(1).tend;
    
end