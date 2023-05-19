% function merge_results reduces results cell array in a Results structure
% with n fields as EigenValue Percent (FFproc Settings). Each component
% will contain m values as Time Windows Evaluated (FFproc Settings).
% version 1.0 / 22apr2020 / cc

function RAW = FFproc_merge_results(results)

    RAW = struct;
    % Concatenating results structure
    % dimensions: (Ntargfreq x %data(eigperc))
    Results = cat(1,results{:});   
    
    for eg = 1:size(Results,2)      % EigenValue Percent

        STR = cat(1,Results(:,eg)); 

        for f = 1:size(Results,1)   % Target Frequencies

            RAW(eg).Zxx = cat(1,STR.Zxx);
            RAW(eg).Zxy = cat(1,STR.Zxy);
            RAW(eg).Zyx = cat(1,STR.Zyx);
            RAW(eg).Zyy = cat(1,STR.Zyy);
            RAW(eg).Zxx_Err = cat(1,STR.Zxx_Err);
            RAW(eg).Zxy_Err = cat(1,STR.Zxy_Err);
            RAW(eg).Zyx_Err = cat(1,STR.Zyx_Err);
            RAW(eg).Zyy_Err = cat(1,STR.Zyy_Err);
            RAW(eg).txz = cat(1,STR.txz);
            RAW(eg).tyz = cat(1,STR.tyz);
            RAW(eg).txz_Err = cat(1,STR.txz_Err);
            RAW(eg).tyz_Err = cat(1,STR.tyz_Err);
            RAW(eg).phi11 = cat(1,STR.phi11);
            RAW(eg).phi12 = cat(1,STR.phi12);
            RAW(eg).phi21 = cat(1,STR.phi21);
            RAW(eg).phi22 = cat(1,STR.phi22);
            RAW(eg).phimin = cat(1,STR.phimin);
            RAW(eg).phimax = cat(1,STR.phimax);
            RAW(eg).alpha = cat(1,STR.alpha);
            RAW(eg).beta = cat(1,STR.beta);
            RAW(eg).theta = cat(1,STR.theta);
            RAW(eg).lambda = cat(1,STR.lambda);
            RAW(eg).phimin_Err = cat(1,STR.phimin_Err);
            RAW(eg).phimax_Err = cat(1,STR.phimax_Err);
            RAW(eg).alpha_Err = cat(1,STR.alpha_Err);
            RAW(eg).beta_Err = cat(1,STR.beta_Err);
            RAW(eg).theta_Err = cat(1,STR.theta_Err);
            RAW(eg).lambda_Err = cat(1,STR.lambda_Err);

        end
    end
            
end