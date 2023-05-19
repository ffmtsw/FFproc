function [Zxx,Zxy,Zyx,Zyy,txz,tyz] = bivariate(N,nts,data,act_chan,magref)

    % Number of EM fields
    nchan = 5;

    % Indices for channels (ex, ey, bx, by, bz)
    ex = 1; ey = 2; bx = 3; by = 4; bz = 5;

    % Preallocating space for EM fields
    EX = zeros(N,nts);      EY = zeros(N,nts);
    BX = zeros(N,nts);      BY = zeros(N,nts);      BZ = zeros(N,nts);
    BXR = zeros(N,nts);     BYR = zeros(N,nts);

    % Pre-allocating space for Z Tensor
    Zxx = NaN(1,nts);       Zxy = NaN(1,nts);
    Zyx = NaN(1,nts);       Zyy = NaN(1,nts);

    % Pre-allocating space for Tipper Vector
    txz = NaN(1,nts);       tyz = NaN(1,nts);

    % Index for active channels
    ind_chan = cat(1,act_chan);
    
    for a = 1:nts

        % Channel jump
        chjump = nchan*(a - 1);
        rchjumpx = nchan*(magref(a,1)-1);
        rchjumpy = nchan*(magref(a,2)-1);

        if ind_chan{a}(ex) && ind_chan{a}(ey)
            % Retrieving horizontal EM-Field components
            EX(:,a) = data(:,ex + chjump);         EY(:,a) = data(:,ey + chjump);
            BX(:,a) = data(:,bx + chjump);         BY(:,a) = data(:,by + chjump);
            % Accounting for remote reference
            BXR(:,a) = data(:,bx + rchjumpx);   
            BYR(:,a) = data(:,by + rchjumpy);
            % Impedance bivariate estimation
            Zxx(a) = bivarZ('XX',EX(:,a),EY(:,a),BX(:,a),BY(:,a),BXR(:,a),BYR(:,a));
            Zxy(a) = bivarZ('XY',EX(:,a),EY(:,a),BX(:,a),BY(:,a),BXR(:,a),BYR(:,a));
            Zyx(a) = bivarZ('YX',EX(:,a),EY(:,a),BX(:,a),BY(:,a),BXR(:,a),BYR(:,a));
            Zyy(a) = bivarZ('YY',EX(:,a),EY(:,a),BX(:,a),BY(:,a),BXR(:,a),BYR(:,a));
        else
            Zxx(a) = complex(NaN,NaN);
            Zxy(a) = complex(NaN,NaN);
            Zyx(a) = complex(NaN,NaN);
            Zyy(a) = complex(NaN,NaN);
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
            txz(a) = bivarT('XZ',BX(:,a),BY(:,a),BZ(:,a),BXR(:,a),BYR(:,a));
            tyz(a) = bivarT('YZ',BX(:,a),BY(:,a),BZ(:,a),BXR(:,a),BYR(:,a));
        else
            % Tipper bivariate estimation
            txz(a) = complex(NaN,NaN);
            tyz(a) = complex(NaN,NaN);
        end
    end

end