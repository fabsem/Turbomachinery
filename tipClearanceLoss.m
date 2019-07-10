function [deltaEta] = tipClearanceLoss(mdot,deltaC,b,c,U,va,leul,beta1,beta2,eta,S,rho)

%TIPCLEARANCELOSS Estimate Tip Clearance losses evaluating the derating in efficiency
%
% INPUT DATA:
% deltaC -->  Clearance height
% b --> Blade height
% c --> Chord
% leul --> Real Work
% eta --> Stage efficiency
%
% Example:
%
%     [deltaEta] = tipClearanceLoss_deltaP(deltaC,b,c,U,leul,eta)
%
% by Lakshiminarayana, 1971

% Average angle
beta_av= (beta1 + beta2) / 2;

% Flow Coefficient
phi= abs(va / U); %pi*mdot / (S*U*rho);

% Work Coefficient
psi= 40;%leul*eta / U^2;

% Derate in efficiency
deltaEta = (0.07 * deltaC / b * psi ) /cosd(beta_av) * ...
           (1 + 10 * sqrt((phi * deltaC / c) / (psi * cosd(beta_av)) ) );

end
