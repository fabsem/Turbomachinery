function [b_lieblein] = lieblein_bcoeff(beta1)

% LIEBLEIN_BCOEFF Returns the b coefficient for deviation angle calculation
%
% Example:
%
%   [b_lieblein] = lieblein_bcoeff(beta1)

    beta1 = abs(beta1);
    
    b_lieblein = 0.9625 - 0.17 * beta1/100 - 0.85 * (beta1/100)^3;

end