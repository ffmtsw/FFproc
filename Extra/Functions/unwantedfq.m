% function remfq removes single frequencies and frequency multiples
% introduced on job parameters.

% version 1.0 / 01jun2014 / ph
% version 2.0 / 24may2019 / at
% version 3.0 / 09may2020 / cc      Modified to match the newest version of
%                                   FFMT 

function remfq = unwantedfq(job,freq,ntseg)

    sing = job.singf;
    mult = job.mult;
    sing_w = job.phsing/2;
    mult_w = job.phmult/2;

    tf = true(size(freq));
    
    if ~isempty(mult)    
        % Multiples removal
        for f = 1:numel(mult)
            m = ceil(max(freq)/mult(f));
            for k = 1:m
                tf(freq >= (k*mult(f))-mult_w & freq < (k*mult(f))+mult_w) = false;
            end
        end
    end
    
    if ~isempty(sing)
        % Single frequencies removal
        for f = 1:numel(sing)
            tf(freq >= sing(f)-sing_w & freq < sing(f)+sing_w) = false;
        end
    end   

    remfq = repmat(tf,ntseg,1)';
    
end