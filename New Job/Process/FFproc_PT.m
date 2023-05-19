% function FFproc_PT calculates the elements of the phase tensor (PT) as
% phi_max, phi_min, alpha, beta, theta and lambda from Impedance Matrix Z

% version 1.1 / 31mar2011 / aj
% version 1.2 / 15jun2011 / aj
% version 2.0 / 22apr2020 / cc  The element theta and lambda have been
%                               added to the output results.

function PT = FFproc_PT(Z)

    X = real(Z);    X = reshape(X,2,2); 
    Y = imag(Z);    Y = reshape(Y,2,2);
    PHI = X\Y;
    
    Pdet = det(PHI);
    PHI1 = trace(PHI)/2;
    PHI3 = (PHI(1,2) - PHI(2,1))/2;
    phimin = abs(sqrt(PHI1^2 + PHI3^2) - sqrt(abs(PHI1^2 + PHI3^2 - Pdet)));
    phimax = sqrt(PHI1^2 + PHI3^2) + sqrt(abs(PHI1^2 + PHI3^2 - Pdet));
    alpha = 0.5*atan2((PHI(1,2)+PHI(2,1)),(PHI(1,1)-PHI(2,2)));
    beta = 0.5*atan2(PHI3,PHI1);
    theta = alpha - beta;
    lambda = (phimax - phimin)./(phimax + phimin);
    
    PT.phi11 = PHI(1,1);    PT.phi12 = PHI(1,2);
    PT.phi21 = PHI(2,1);    PT.phi22 = PHI(2,2);
    PT.alpha = alpha;       PT.beta = beta;     
    PT.theta = theta;       PT.lambda = lambda;
    PT.phimax = phimax;     PT.phimin = phimin;

end