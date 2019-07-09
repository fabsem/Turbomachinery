function [Y] = profileLossesLieblein(theta_c,sigma,beta1,beta2,varargin)

%PROFILELOSSESLIEBLEIN Estimate profile losses using Lieblein correlation.
%
% Beware: It has to be treated adding Re, Ma, thickness corrections.

    if ~isempty(varargin)
        H = varargin;

    else

        H = 1.08;


    end

    Y = (2 + theta_c * sigma / cosd(beta2) * (cosd(beta1)/cosd(beta2))^2 * (2 * H / (3 * H - 1)))/...
            (1 - theta_c * sigma / cosd(beta2) * H)^3;
end
