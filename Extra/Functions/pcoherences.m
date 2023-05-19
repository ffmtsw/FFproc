function [rZxx,rZxy,rZyx,rZyy,rtxz,rtyz] = pcoherences(N,nts,data,act_chan,magref)

    % Number of EM fields
    nchan = 5;

    % Indices for channels (ex, ey, bx, by, bz)
    ex = 1; ey = 2; bx = 3; by = 4; bz = 5;

    % Preallocating space for EM fields
    EX = zeros(N,nts);    EY = zeros(N,nts);
    BX = zeros(N,nts);    BY = zeros(N,nts);    BZ = zeros(N,nts);
    BXR = zeros(N,nts);   BYR = zeros(N,nts);

    % Pre-allocating space for Z Tensor
    rZxx = zeros(1,nts);      rZxy = zeros(1,nts);
    rZyx = zeros(1,nts);      rZyy = zeros(1,nts);

    % Pre-allocating space for Tipper Vector
    rtxz = zeros(1,nts);      rtyz = zeros(1,nts);

    % Index for active channels
    ind_chan = cat(1,act_chan);

    for a = 1:nts       

        % Channel jump
        chjump = nchan*(a - 1);
        rchjumpx = nchan*(magref(a,1)-1);
        rchjumpy = nchan*(magref(a,2)-1);

        if ind_chan{a}(ex) && ind_chan{a}(ey)
            % Retrieving horizontal EM-Field components
            EX(:,a) = data(:,ex + chjump);    EY(:,a) = data(:,ey + chjump);
            BX(:,a) = data(:,bx + chjump);    BY(:,a) = data(:,by + chjump);
            % Accounting for remote reference
            BXR(:,a) = data(:,bx + rchjumpx);   
            BYR(:,a) = data(:,by + rchjumpy);
            % Impedance bivariate estimation
            % Calculation of partial coherences
            rZxx(a) = abs(partialcorr(EX(:,a),BX(:,a),BY(:,a)));
            if isnan(rZxx(a))
                rZxx(a) = 0;
            end
            rZxy(a) = abs(partialcorr(EX(:,a),BY(:,a),BX(:,a)));
            if isnan(rZxy(a))
                rZxy(a) = 0;
            end
            rZyx(a) = abs(partialcorr(EY(:,a),BX(:,a),BY(:,a)));
            if isnan(rZyx(a))
                rZyx(a) = 0;
            end
            rZyy(a) = abs(partialcorr(EY(:,a),BY(:,a),BX(:,a)));
            if isnan(rZyy(a))
                rZyy(a) = 0;
            end
        else
            rZxx(a) = 0;
            rZxy(a) = 0;
            rZyx(a) = 0;
            rZyy(a) = 0;
        end
        
        if ind_chan{a}(bz)
            if all(isnan(BX(:,a))) || all(isnan(BY(:,a)))
                % Retrieving horizontal magnetic fields
                BX(:,a) = data(:,bx + chjump);    BY(:,a) = data(:,by + chjump);                
                % Accounting for remote reference
                BXR(:,a) = data(:,bx + rchjumpx);   
                BYR(:,a) = data(:,by + rchjumpy);      
            end
            % Retrieving vertical magnetic fields
            BZ(:,a) = data(:,bz + chjump);
            % Tipper bivariate estimation
            rtxz(a) = abs(partialcorr(BZ(:,a),BX(:,a),BY(:,a)));
            if isnan(rtxz(a))
                rtxz(a) = 0;
            end
            rtyz(a) = abs(partialcorr(BZ(:,a),BY(:,a),BX(:,a)));
            if isnan(rtyz(a))
                rtyz(a) = 0;
            end
        else
            % Tipper bivariate estimation
            rtxz(a) = 0;
            rtyz(a) = 0;
        end
    end

end