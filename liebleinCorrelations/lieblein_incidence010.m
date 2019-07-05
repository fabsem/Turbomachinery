function [i010] = lieblein_incidence010(beta1,sigma)

% LIEBLEIN_INCIDENCE010 Returns the reference incidence for symmetric
% profile, referred to the max thickness case of 10% of the chord
%
% Example:
%
%   [i010] = lieblein_incidence010(beta1,sigma)

    beta1 = abs(beta1);
    
    i010 = (beta1^(0.914 + (sigma^3)/160))/(5 + 46 * exp(-2.3 * sigma)) - ...
            0.1 * sigma^3 * exp((beta1 - 70)/4);

end