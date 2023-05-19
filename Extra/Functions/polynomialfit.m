function [comp,comp_err,id_re,id_im,cf_re,cf_im] = polynomialfit(freq,arg,err,argid,polyord)

    % Logical flag to use either TFcomponent or log10(TFcomponent)
    tflog = false;
    
    % Period vector
    per = 1./freq;

    switch upper(argid)
        case {'ZXY','ZYX'}
            % Real Part
            s = nanmedian(sign(real(arg(:))));
            [re,err_re,id_re,cf_re] = fitpoly2(log10(per),log10(abs(real(arg))),err,polyord);
            comp_re = s.*10.^re;        
            % Imaginary Part
            s = nanmedian(sign(imag(arg(:))));
            [im,err_im,id_im,cf_im] = fitpoly2(log10(per),log10(abs(imag(arg))),err,polyord);
            comp_im = s.*10.^im;        
            % Complex component
            comp = complex(comp_re(:),comp_im(:));
            comp(isnan(real(comp))|isnan(imag(comp))) = complex(NaN,NaN);        
            % Error
            comp_err = mean([err_re(:),err_im(:)],2);    
            comp_err(isnan(real(comp))|isnan(imag(comp))) = NaN;
        case {'ZXX','ZYY'}
            if tflog
                % Real Part
                s = nanmedian(sign(real(arg(:)))); %#ok<UNRCH,*NANMEDIAN> 
                [re,err_re,id_re,cf_re] = fitpoly2(log10(per),log10(abs(real(arg))),err,polyord);
                comp_re = s.*10.^re;        
                % Imaginary Part
                s = nanmedian(sign(imag(arg(:))));
                [im,err_im,id_im,cf_im] = fitpoly2(log10(per),log10(abs(imag(arg))),err,polyord);
                comp_im = s.*10.^im;  
            else
                % Real Part
                [comp_re,err_re,id_re,cf_re] = fitpoly2(log10(per),real(arg),err,polyord);
                % Imaginary Part
                [comp_im,err_im,id_im,cf_im] = fitpoly2(log10(per),imag(arg),err,polyord); 
            end
            % Complex component
            comp = complex(comp_re(:),comp_im(:));
            comp(isnan(real(comp))|isnan(imag(comp))) = complex(NaN,NaN); 
            % Error
            comp_err = mean([err_re(:),err_im(:)],2);    
            comp_err(isnan(real(comp))|isnan(imag(comp))) = NaN;           
        case {'TXZ','TYZ'}
            if all(isnan(arg(:)))
                comp = NaN(numel(freq),1);            
                comp_err = NaN(numel(freq),1);
            else
                [comp_re,err_re,id_re,cf_re] = fitpoly2(log10(per),real(arg),err,polyord);
                [comp_im,err_im,id_im,cf_im] = fitpoly2(log10(per),imag(arg),err,polyord);
                % Complex component
                comp = complex(comp_re(:),comp_im(:));
                comp(isnan(real(comp))|isnan(imag(comp))) = complex(NaN,NaN);
                % Error
                comp_err = mean([err_re(:),err_im(:)],2);    
                comp_err(isnan(real(comp))|isnan(imag(comp))) = NaN; 
            end
    end

end