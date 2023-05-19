% shrink takes a distribution and re-define weights, from 0 to 1 based on
% the ammount of data it needs to be considered.
% It is assumed that data is from [0 , 1].

function w = shrink(w,rnk,perc)

    % Number of data
    N = numel(w);

    % Ammount of data to be considered
    maxw = ceil(rnk*N/100);

    % If only one datum can be used, assign 1 , otherwise assign 0
    if maxw == 1
        w = zeros(size(w));
        w(1) = 1;
    elseif numel(unique(w)) == 1
        w = ones(size(w));    
    else        
        % Data outbounds should be tagged as zero
        w(maxw + 1:end) = 0;
        % Find unique vector
        wunique = unique(w,'stable');
        % Find the minimum value (will be the new zero)
        if numel(wunique) <= 2
           minw = wunique(end);
        else
             minw = wunique(end-1);
        end
        % Tag all the zero values to be the minimum
        w(w==0) = minw;
        % Shift distribution to make the minimum -> zero
        w = w-minw;
        % Distribution normalization (from 0 to 1)
        w = w./max(w);
    end

    if perc > 0 && N > 2
        n = ceil(N*perc/100);
        w(1:n) = 0;
        w(w<0) = 0;
        w(isnan(w) | isinf(w)) = 0;
        w = shrink(w,rnk,100);
    end

end