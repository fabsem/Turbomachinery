function [n_lieblein] = lieblein_ncoeff(beta1,sigma)

% LIEBLEIN_NCOEFF Returns the n coefficient for incidence angle calculation
%
% Example:
%
%   [Ki_thick] = lieblein_icorrectionThick(t,c)


    beta1 = abs(beta1);
    
    n_lieblein = 0.025 * sigma - 0.06 - ((beta1/90)^(1 + 1.2 * sigma))/(1.5 + 0.43 * sigma);

end