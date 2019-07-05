function [delta010] = lieblein_delta010(beta1,sigma)

% LIEBLEIN_DELTA010 Returns the reference incidence for symmetric
% profile, referred to the max thickness case of 10% of the chord
%
% Example:
%
%   [delta010] = lieblein_delta010(beta1,sigma)

    beta1 = abs(beta1);
    
    delta010 = 0.01 * sigma * beta1 + (0.74 * sigma ^ 1.9 + 3 * sigma) * (beta1/90) ^ (1.67 + 1.09 * sigma);

end