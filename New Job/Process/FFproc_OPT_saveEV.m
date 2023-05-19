function EV = FFproc_OPT_saveEV(job,tfreq_vect,n_tseg,t_seg,id,EV)  

    for f = 1:numel(tfreq_vect)
        EV(f).perc = job.ranking;
        EV(f).sr = job.sr;
        EV(f).tfreq = tfreq_vect;
        EV(f).n_tseg = n_tseg;
        EV(f).tseg = t_seg;
        EV(f).tstart = job.ffts(1).tstart;
        EV(f).tend = job.ffts(1).tend;
        if numel(fieldnames(id)) == 0
            for a = 1:numel(job.ffts)
                EV(f).Zxx_OPT(a).id = [];
                EV(f).Zxy_OPT(a).id = [];
                EV(f).Zyx_OPT(a).id = [];
                EV(f).Zyy_OPT(a).id = [];
                EV(f).txz_OPT(a).id = [];
                EV(f).tyz_OPT(a).id = [];                
            end
        else
            for a = 1:numel(job.ffts)
                EV(f).Zxx_OPT(a).id = [id(a).Zxx_Re(:),id(a).Zxx_Im(:)];
                EV(f).Zxy_OPT(a).id = [id(a).Zxy_Re(:),id(a).Zxy_Im(:)];
                EV(f).Zyx_OPT(a).id = [id(a).Zyx_Re(:),id(a).Zyx_Im(:)];
                EV(f).Zyy_OPT(a).id = [id(a).Zyy_Re(:),id(a).Zyy_Im(:)];
                EV(f).txz_OPT(a).id = [id(a).txz_Re(:),id(a).txz_Im(:)];
                EV(f).tyz_OPT(a).id = [id(a).tyz_Re(:),id(a).tyz_Im(:)];                
            end
        end
    end
    
end