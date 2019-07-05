function [Ki_thick] = lieblein_icorrectionThick(t_c)

% LIEBLEIN_ICORRECTIONTHICK Returns the thickness parameter for incidence
% angle correction
%
% Example:
%
%   [Ki_thick] = lieblein_icorrectionThick(t_c)


    Ki_thick = (10 * t_c)^(0.28 /(0.1 + (t_c)^0.3));

end