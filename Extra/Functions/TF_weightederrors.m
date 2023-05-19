function [Zwmed,Z_Err] = TF_weightederrors(Z,N,W)

    % Weighted median for transfer function
    Zwmed = complex(wmedian(real(Z),W),wmedian(imag(Z),W));
    % Errors after Häuserer(2010)
    dZ = max([1.483*wmedian(abs(real(Z) - real(repmat(Zwmed,N,1))),W),...
              1.483*wmedian(abs(imag(Z) - imag(repmat(Zwmed,N,1))),W)]);    
    Z_Err = 1.96*dZ./sqrt(N);

end