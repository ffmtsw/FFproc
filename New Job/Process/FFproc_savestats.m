function stat = FFproc_savestats(n_tseg,n_samples,t_seg,osc,n_skip,overlap_tf,time)
    
    stat.n_tseg = n_tseg;
    stat.n_samples = n_samples;
    stat.t_seg = t_seg;
    stat.osc = osc;
    stat.n_skip = n_skip;
    stat.overlap_tf = overlap_tf;
    stat.time = time;

end