function [M2,P2,T2,Pt2,rho2,Yshock]=shockLosses(M1,Pt1,P1,T1,rho1)
gamma=1.4;

P2 = P1 * ( 2*gamma*M1^2-(gamma-1) ) / (gamma+1);
T2 = T1 * (2*gamma*M1^2-gamma+1)*((gamma-1)*M1^2+2)/((gamma+1)^2*M1^2);
rho2 = rho1 * (gamma+1)*M1^2/((gamma-1)*M1^2+2);
Pt2 = Pt1 * ((gamma+1)*M1^2/((gamma-1)*M1^2+2) )^(gamma/(gamma-1)) * ((gamma+1)/(2*gamma*M1^2-gamma+1))^(1/(gamma-1));
%Tt2 = Tt1;
M2 = ((gamma-1)*M1^2+2)/(2*gamma*M1^2-gamma+1);

Yshock=(Pt1-Pt2)/(Pt1-P1);

end