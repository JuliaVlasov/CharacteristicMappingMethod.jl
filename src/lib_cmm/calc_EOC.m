function [EOC] = calc_EOC(u_ref, u_h1, u_h2, h1, h2, ord, weights)
% This function calculates the experimental order of convergence.
% See for explanation: https://de.wikipedia.org/wiki/Experimentelle_Konvergenzordnung

    if (nargin <= 5) || isempty(ord)
        ord = "inf";
    end
    if (nargin <= 6) || isempty(weights)
        weights = 1;
    end
    
    if isempty(u_h1)
        EOC = log10(norm(u_ref,ord)/norm(u_h2,ord))/log10(h1/h2);
    else
        EOC = log10(norm((u_h1-u_ref).*weights,ord)/norm((u_h2-u_ref).*weights,ord))/log10(h1/h2);
    end
end