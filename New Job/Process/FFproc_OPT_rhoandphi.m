function FFsave = FFproc_OPT_rhoandphi(FFsave,job)

    for i = 1:numel(FFsave)        
     
        % Resistivity, Phase and Errors
        [rho,phi] = rhoandphi(FFsave(i).Zxx,FFsave(i).Zxy,...
                              FFsave(i).Zyx,FFsave(i).Zyy,...
                              FFsave(i).Zxx_Err,FFsave(i).Zxy_Err,...
                              FFsave(i).Zyx_Err,FFsave(i).Zyy_Err,...
                              FFsave(i).freq);
                          
        % Transfer Function scaling and correction
        % Resistivity
        FFsave(i).rhoxx = rho.rhoxx;
        FFsave(i).rhoxy = rho.rhoxy;
        FFsave(i).rhoyx = rho.rhoyx;
        FFsave(i).rhoyy = rho.rhoyy;
         % Resistivity Errors
        FFsave(i).rhoxx_Err = rho.rhoxx_Err;
        FFsave(i).rhoxy_Err = rho.rhoxy_Err;
        FFsave(i).rhoyx_Err = rho.rhoyx_Err;
        FFsave(i).rhoyy_Err = rho.rhoyy_Err;
        % Phase
        FFsave(i).phixx = phi.phixx;
        FFsave(i).phixy = phi.phixy;                
        FFsave(i).phiyx = phi.phiyx;        
        FFsave(i).phiyy = phi.phiyy;
        % Phase Errors
        FFsave(i).phixx_Err = phi.phixx_Err;
        FFsave(i).phixy_Err = phi.phixy_Err;
        FFsave(i).phiyx_Err = phi.phiyx_Err;
        FFsave(i).phiyy_Err = phi.phiyy_Err;  

    end
    
end