% wmedian function 
% x input: data matrix. The wmedian will be calculated for each column of x
% 
% w input: weighted function contains weights [0,1]. 
%          The weights are considered as the times its value will be taken
%          into account for median calculation:
% W = 1 means 100 times
% W = 0.5 means 50 times
% W = 0.1 means 10 times
% The median is calculated only for those values different from NaN
%
% version 1.0 / 16aug2022 / cc

function xmed = wmedian(x,w)
    
    try
        [sortx,order] = sort(x,1);
        sortw = w(order);
    
        midpoint = sum(sortw)/2;
        csumw = cumsum(sortw,1);
        
        for i = 1:length(midpoint)
            j = find(csumw(:,i) <= midpoint(i),1,'last');
            if csumw(j,i) == midpoint(i)
                 xmed(i) = mean(sortx([j j+1],i));
            else
                 xmed(i) = sortx(j+1,i);
            end
        end
    catch
        [~,col] = size(x);
        % How many times the value should be taken into account for median 
        % estimation
        W = round(w*100);
        % Reducing weights less than zero;
        W(W<0) = 0;

        xmed = NaN(col,1);
        for i = 1:col
            % Repeat W times the element
            X = repelem(x(:,i),W);
            % Calculate median from X (only for non-NaN entries)
            xmed(i) = median(X(~isnan(X)));
        end
    end
    
end