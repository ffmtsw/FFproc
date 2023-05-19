function RAW = FFproc_RAW_rhoandphi(RAW,job)

    for i = 1:numel(RAW)
        for j = 1:numel(job.name)            
            % Resistivity, Phase and Errors
            try
                [rho,phi] = rhoandphi(RAW(i).Zxx(:,j),RAW(i).Zxy(:,j),...
                                      RAW(i).Zyx(:,j),RAW(i).Zyy(:,j),...
                                      RAW(i).Zxx_Err(:,j),RAW(i).Zxy_Err(:,j),...
                                      RAW(i).Zyx_Err(:,j),RAW(i).Zyy_Err(:,j),...
                                      RAW(i).freq);
                                  
                RAW(i).rhoxx(:,j) = rho.rhoxx;
                RAW(i).rhoxy(:,j) = rho.rhoxy;
                RAW(i).rhoyx(:,j) = rho.rhoyx;
                RAW(i).rhoyy(:,j) = rho.rhoyy;
                RAW(i).rhoxx_Err(:,j) = rho.rhoxx_Err;
                RAW(i).rhoxy_Err(:,j) = rho.rhoxy_Err;
                RAW(i).rhoyx_Err(:,j) = rho.rhoyx_Err;
                RAW(i).rhoyy_Err(:,j) = rho.rhoyy_Err;
                RAW(i).phixx(:,j) = phi.phixx;
                RAW(i).phixy(:,j) = phi.phixy;
                RAW(i).phiyx(:,j) = phi.phiyx;
                RAW(i).phiyy(:,j) = phi.phiyy;
                RAW(i).phixx_Err(:,j) = phi.phixx_Err;
                RAW(i).phixy_Err(:,j) = phi.phixy_Err;
                RAW(i).phiyx_Err(:,j) = phi.phiyx_Err;
                RAW(i).phiyy_Err(:,j) = phi.phiyy_Err;
            catch
                RAW(i).rhoxx(:,j) = NaN;
                RAW(i).rhoxy(:,j) = NaN;
                RAW(i).rhoyx(:,j) = NaN;
                RAW(i).rhoyy(:,j) = NaN;
                RAW(i).rhoxx_Err(:,j) = NaN;
                RAW(i).rhoxy_Err(:,j) = NaN;
                RAW(i).rhoyx_Err(:,j) = NaN;
                RAW(i).rhoyy_Err(:,j) = NaN;
                RAW(i).phixx(:,j) = NaN;
                RAW(i).phixy(:,j) = NaN;
                RAW(i).phiyx(:,j) = NaN;
                RAW(i).phiyy(:,j) = NaN;
                RAW(i).phixx_Err(:,j) = NaN;
                RAW(i).phixy_Err(:,j) = NaN;
                RAW(i).phiyx_Err(:,j) = NaN;
                RAW(i).phiyy_Err(:,j) = NaN;
            end            
        end
    end

end