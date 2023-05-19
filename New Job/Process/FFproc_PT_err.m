% function FFproc_PT_err calculates the error from the elements of the 
% phase tensor (PT): phi_max_sig, phi_min_sig, alpha_sig, beta_sig, 
% theta_sig and lambda_dig from Impedance Matrix Z

% version 1.1 / 31mar2011 / aj
% version 2.0 / 22apr2020 / cc  The errors calculations for the elements 
%                               theta and lambda have been added to the 
%                               output results.

function PT_err = FFproc_PT_err(Z,Z_err)

    N = 1000;
    PTr = zeros(6,N);
    for i = 1:N
       Zr =  Z + complex(randn(4,1),randn(4,1)).*Z_err;
       PT = FFproc_PT(Zr);
       PTr(1:2,i) = mod([PT.phimin,PT.phimax],pi/2);
       PTr(3:4,i) = mod([PT.beta,PT.alpha],pi);
       PTr(5,i) = mod(PTr(4,1) - PTr(3,1),pi);
       PTr(6,i) = mod((PT.phimax - PT.phimin)./(PT.phimax + PT.phimin),pi);    
    end
    
    med = median(PTr,2)';
    PT_err = median(abs(PTr'-ones(N,1)*med))*1.4;

end