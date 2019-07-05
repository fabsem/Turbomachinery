function [Kdelta_thick] = lieblein_deltacorrectionThick(t_c)

% LIEBLEIN_ICORRECTIONTHICK Returns the thickness parameter for deviation
% angle correction
%
% Example:
%
%   [Kdelta_thick] = lieblein_deltacorrectionThick(t_c)


    Kdelta_thick = 6.25 * t_c + 37.5 * t_c^2;
    
end