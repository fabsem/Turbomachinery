function [deltaPt] = tipClearanceLoss_deltaP(Z,Nrow,Gamma,rho1,rho2,deltaC,r1,r2,Vm1,Vm2,Vt1,Vt2,mdot)

%TIPCLEARANCELOSS Estimate Tip Clearance losses evaluating the torque and the leakage mass flow
%
% INPUT DATA:
% (1 = entrance, 2 = exit)
% Z --> Number of blades in the row
% deltaC -->  Clearance height
% r1 --> Tip radius
% rho --> Flow Density
% Vm --> Meridional (axial) velocity
% Vt --> Tangential velocity
% c --> Chord
% Gamma --> Stagger angle (beta_in+i+theta/2)
% rho_average --> mean(rho1,rho2)
% tau ---> Torque due to clearance flow
% Nrow --> Number of the row, first rotor or second
%
% Example:
%
%     [zetaTC] = tipClearanceLoss()
%
% by Aungier - Axial Flow Compressors (2003) p. 146

rho_average = (rho1 + rho2) / 2;

% Blade Torque due to clearance flow
tau = abs(pi * deltaC * ((r1 * rho1 * Vm1) + (r2 * rho2 * Vm2)) * (r2 * Vt2 - r1 * Vt1));

% Average Pressure Difference across each blade in the blade row
deltaP = tau / (Z * (r1 + r2)/2 * deltaC * cosd(Gamma));

% Velocity of the leakage flow
Uc = 0.816 * sqrt(2 * deltaP / rho_average) / Nrow^0.2;

% Leakage Mass flow rate
mdotC = rho_average * Uc * Z * deltaC * cosd(Gamma);

%Total pressure loss for the entire blade row
deltaPt = deltaP * mdotC / mdot;


end
