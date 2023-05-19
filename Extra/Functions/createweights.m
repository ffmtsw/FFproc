function W = createweights(N,n,win,perc,w0,pow)

    % N: Number of all elements (ntseg)
    % n: Number of selected elements (seg_eval)    
    % l: Width of window (only for chi function)
    
    if isempty(win)
        win = 'cos';
        pow = 1;
    end

    switch win
        case 'flat'
            x = zeros(N,1);
            x(1:n,1) = 1;
            W = x;
            
        case 'ramp'
            x = zeros(N,1);
            x(1:n,1) = linspace(1,0,n);
            W = x;
            
        case 'cos'
            x = zeros(N,1);
            c = tukeywin(n*2,1).^pow;
            x(1:n,1) = c(end/2 + 1:end);
            W = x;

        case 'ham'
            x = zeros(N,1);
            c = hamming(n*2).^pow;
            x = x + c(end);
            x(1:n,1) = c(end/2 + 1:end);
            W = x;

        case 'chi'
            % Length of weighting function
            x = zeros(N,1);
            % m: Number of donweighted elements at the beginning
            m = ceil(perc*n/100);
            % l: Length of chi function
            l = n - m;
            % First part of the chi window
            x1 = (1-w0)*(tukeywin(m*2,1).^1) + w0;
            x1 = x1(1:end/2);
            % Second part of the chi window
            x2 = tukeywin(l*2,1).^2;
            x2 = x2(end/2 + 1:end);

            x(1:m,1) = x1;
            x(m+1:n) = x2;
            W = x;
    end

end