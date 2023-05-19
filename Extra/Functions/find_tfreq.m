% function target_freq results in a logaritmic (centered) equal-spaced 
% vector from the lowest to the highest frequency typed in 
% Fproc_Parameters GUI
% version 1.0 / 08abr2020 / cc

function [tfreq_lim,tfreq_vect,tfreq_cutoff] = find_tfreq(nfreq,maxfreq,minfreq)
    
    % Maximum and minimum log10 of frequencies
    maxfreqlog = ceil(log10(maxfreq));
    minfreqlog = floor(log10(minfreq));

    % Equal spaced log10 frequency vector
    fvect = sort(minfreqlog-1:1/nfreq:maxfreqlog+1,'descend')';
    
    tfreq_lim = [fvect(1:end-1),fvect(2:end)];
    tfreq_vect = mean(tfreq_lim,2);
    tfreq_cutoff = [tfreq_vect(1:end-2),tfreq_vect(3:end)];

    tfreq_lim = 10.^tfreq_lim;
    tfreq_vect = 10.^tfreq_vect;
    tfreq_cutoff = 10.^tfreq_cutoff;

    % Deleting data out of bounds [minfreq,maxfreq]    
    del_ind = (tfreq_vect > maxfreq | tfreq_vect < minfreq);
    tfreq_vect(del_ind,:) = [];
    tfreq_lim(del_ind,:) = [];
    tfreq_cutoff(del_ind(2:end-1),:) = [];
    
end