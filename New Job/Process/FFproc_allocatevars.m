function [n_tseg,n_samples,t_seg,osc,n_skip,overlap_tf,...
          linreg,time,FC,TF,RNK,EigVal,results] = FFproc_allocatevars(tfreq)

    % Preallocating variables (vector)
    n_tseg = NaN(numel(tfreq),1);
    n_samples = NaN(numel(tfreq),1);
    t_seg = NaN(numel(tfreq),1);
    osc = NaN(numel(tfreq),1);
    n_skip = NaN(numel(tfreq),1);
    overlap_tf = NaN(numel(tfreq),1);
    % Preallocating variables (cell)    
    linreg = cell(1,numel(tfreq));
    time = cell(1,numel(tfreq));
    FC = cell(1,numel(tfreq));
    TF = cell(1,numel(tfreq));
    RNK = cell(1,numel(tfreq));
    results = cell(1,numel(tfreq));
    EigVal = cell(1,numel(tfreq));

end