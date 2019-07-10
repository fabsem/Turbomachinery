function [deltaEta] = deltaEta_calc(variable,Y,zeta,M2,P1,P2,T1,T2,Beta,Pt1,v1)

% DELTAETA_CALC Estimate derate in efficiency from Y or Zeta formulations
%
%

gamma = 1.4;
cp = 1004.69;

h1 = cp * T1;
h2is = cp * T1 * Beta^((gamma-1) / gamma);
Pt2 = Pt1 - (Pt1 - P1) * Y;


if strcmp(variable,'Y')

  zeta = (1 - (1 + Y * (1 - P2/Pt2))^((1 - gamma)/gamma)) / ((gamma - 1) / 2 * M2 ^ 2);

end

  h2 = v1^2 / 2 * zeta + h2is;
  etaSTAR = (h2is - h1) / (h2 - h1);
  deltaEta = 1 - etaSTAR;

end
