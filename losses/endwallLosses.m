function [zetaEW] = endwallLosses(zetaP,p2,pt2,M2,deltaBeta,rho1,rho2,V1,V2,alpha2,t_TE,b,S)

% ENDWALLLOSSES Estimate Endwall Losses by means of Lakshiminarayana Correlation
%
% zetaP --> Profile Losses
% epsilon --> Turning angle [rad]
% rho --> Flow density
% V --> Flow Velocity
% S --> Flow Cross Surface
% alpha2 --> Flow Angle
% t_TE --> Trailing Edge Thickness
% b --> blade height
%
% Lakshiminarayana, Fluid Dynamics and Heat Transfer of Turbomachinery (1996) - p. 581

gamma = 1.4;

%zetaP = (1 - (1 + YP * (1 - p2/pt2))^((1 - gamma)/gamma)) / ((gamma - 1) / 2 * M2 ^ 2);

zetaEW = abs(zetaP * (1 + (4 * deltaBeta)/sqrt(rho2 * V2 / (rho1 * V1))) * (S * cosd(alpha2) - t_TE) / b);

% Scale zetaEW to absolute value from percentage
zetaEW = zetaEW * 0.01;
end
