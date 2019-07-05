function [Kdelta_shape] = lieblein_deltacorrectionShape(shape)

% LIEBLEIN_DELTACORRECTIONSHAPE Returns the shape parameter for deviation angle correction
%
% Known Airfoil Parameterizations:
% 
%   1) naca65
%   2) circulararc
%   3) dca
%
% Example:
%
%   [Kdelta_shape] = lieblein_deltacorrectionShape('naca65')

    switch shape
        
        case 'naca65'
            Kdelta_shape = 1;
        case 'circulararc'
            Kdelta_shape = 1.1;
        case 'dca'
            Kdelta_shape = 0.7;
        otherwise
            Kdelta_shape = 1;
end