function [Ki_shape] = lieblein_icorrectionShape(shape)

% LIEBLEIN_ICORRECTIONSHAPE Returns the shape parameter for incidence angle correction
%
% Known Airfoil Parameterizations:
% 
%   1) naca65
%   2) circulararc
%   3) dca
%
% Example:
%
%   [Ki_shape] = lieblein_icorrectionShape('naca65')

    switch shape
        
        case 'naca65'
            Ki_shape = 1;
        case 'circulararc'
            Ki_shape = 1.1;
        case 'dca'
            Ki_shape = 0.7;
        otherwise
            Ki_shape = 1;
end