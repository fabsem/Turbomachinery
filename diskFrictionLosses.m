function [zetaD] = diskFrictionLosses(delta,Dhub,mdot,S,U1,U2,rho1,b,Cp,HUB,eta1,T1)

%DISKFRICTIONLOSSES Estimate Disk Friction Losses that arise when a surface is rotating with respect to the other
%
% Input Parameters:
%
% Re --> Reynolds number
% delta -->  gap size
% R --> Disk Radius
%

%Cm --> Disk Friction Coefficient
%FAI PLOT
Cm = 1.5;

%d = diameter
%U = peripheralSpeed
%rho = density

w_is = 2 * sqrt(2 * ( - eta1 * Cp * (HUB.T2 - T1) + HUB.w1^2 / 2));
U_max = U1(1) + U2(1);

%Flow Coefficient
PHI = mdot / (S * U1(1) * rho1);

%Work Coefficient
PSI = 2 * w_is / U_max^2;

  zetaD = (8 * Cm * Dhub) / (pi * PHI * PSI * b);

end
