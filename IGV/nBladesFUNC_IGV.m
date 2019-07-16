function [HUB, MID, TIP, nBlades] = nBladesFUNC_IGV(HUB,MID,TIP,rho,mu,Re,solidityBest,Dhub,Dmid,Dtip)

%NBLADESFUNC To calculate the number of blades
%
%
% example:
%
% [HUB, MID, TIP, nBlades] = nBladesFUNC(HUB,MID,TIP,rho1,mu,1e6,solidityBest,Dhub,Dmid,Dtip)

  HUB.c_IGV = Re * mu / (rho * HUB.v1a);
  %HUB.c1 = 0.15;
  MID.c_IGV = Re * mu / (rho * MID.v1a);
  %MID.c1 = 0.10;
  TIP.c_IGV = Re * mu / (rho * TIP.v1a);

  % Blade number
  MID.sigma_IGV = solidityBest(2);
  MID.s_IGV = MID.c_IGV / MID.sigma_IGV;

  nBlades1 = pi * Dmid / MID.s_IGV;
  nBlades1 = floor(nBlades1);
  if rem(nBlades1,2)
      nBlades1 = nBlades1 + 1;
  end

  HUB.s_IGV = pi * Dhub / nBlades1;
  MID.s_IGV = pi * Dmid / nBlades1;
  TIP.s_IGV = pi * Dtip / nBlades1;

  %MID.sigma1 IMPOSED;
  TIP.sigma_IGV =  TIP.c_IGV / TIP.s_IGV;
  HUB.sigma_IGV = HUB.c_IGV / HUB.s_IGV;

  %% Camber angle --> need iteration
  % epsilon = deltaBeta = theta + i - delta
  HUB.epsilon_IGV = HUB.deltaAlfa;
  MID.epsilon_IGV = MID.deltaAlfa;
  TIP.epsilon_IGV = TIP.deltaAlfa;

  %guess
  HUB.theta_IGV = HUB.epsilon_IGV;
  MID.theta_IGV = MID.epsilon_IGV;
  TIP.theta_IGV = TIP.epsilon_IGV;

  nBlades = nBlades1;

  HUB.Re_IGV = rho * HUB.v1a * Dhub / mu;
  MID.Re_IGV = rho * MID.v1a * Dhub / mu;
  TIP.Re_IGV = rho * TIP.v1a * Dhub / mu;


end
