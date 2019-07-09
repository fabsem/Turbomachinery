function [m_lieblein] = lieblein_mcoeff(beta1,shape)

% LIEBLEIN_MCOEFF Returns the n coefficient for deviation angle calculation
%
% Example:
%
%   [m_lieblein] = lieblein_mcoeff(beta1,shape)


    beta1 = abs(beta1);

    switch shape

        case 'naca65'

            m_lieblein = 0.17 - 0.0333 * beta1/100 + 0.333 * (beta1/100)^2;

        case 'DCA'

            m_lieblein = 0.249 + 0.074 * beta1/100 - 0.132 * (beta1/100)^2 + 0.316 * (beta1/100)^3;

    end
end
