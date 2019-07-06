function [HUB, MID, TIP, nBlades] = nBladesFUNC(HUB,MID,TIP,rotor,rho,mu,Re,solidityBest,Dhub,Dmid,Dtip)

%NBLADESFUNC To calculate the number of blades
%
% choose the rotor under analysis.
%
% example:
%
% [HUB, MID, TIP, nBlades] = nBladesFUNC(HUB,MID,TIP,1,rho1,mu,1e6,solidityBest,Dhub,Dmid,Dtip)

if rotor==1
  HUB.c1 = Re * mu / (rho * HUB.w1);
  %HUB.c1 = 0.15;
  MID.c1 = Re * mu / (rho * MID.w1);
  %MID.c1 = 0.10;
  TIP.c1 = Re * mu / (rho * TIP.w1);

  % Blade number
  MID.sigma1 = solidityBest(2);
  MID.s1 = MID.c1 / MID.sigma1;

  nBlades1 = pi * Dmid / MID.s1;
  nBlades1 = floor(nBlades1);
  if rem(nBlades1,2)
      nBlades1 = nBlades1 + 1;
  end

  HUB.s1 = pi * Dhub / nBlades1;
  MID.s1 = pi * Dmid / nBlades1;
  TIP.s1 = pi * Dtip / nBlades1;


  %MID.sigma1 IMPOSED;
  TIP.sigma1 =  TIP.c1 / TIP.s1;
  HUB.sigma1 = HUB.c1 / HUB.s1;

  %% Camber angle --> need iteration
  % epsilon = deltaBeta = theta + i - delta
  HUB.epsilon1 = HUB.deltaBeta1;
  MID.epsilon1 = MID.deltaBeta1;
  TIP.epsilon1 = TIP.deltaBeta1;

  %guess
  HUB.theta1 = HUB.epsilon1;
  MID.theta1 = MID.epsilon1;
  TIP.theta1 = TIP.epsilon1;

  nBlades = nBlades1;

  HUB.Re1 = rho * HUB.w1 * Dhub / mu;
  MID.Re1 = Re;
  TIP.Re1 = Re;

elseif rotor==2
  HUB.c2 = Re * mu / (HUB.rho2 * HUB.w3);
  MID.c2 = Re * mu / (MID.rho2 * MID.w3);
  TIP.c2 = Re * mu / (TIP.rho2 * TIP.w3);

  % Blade number
  MID.sigma2 = solidityBest(5);
  MID.s2 = MID.c2 / MID.sigma2;

  nBlades2 = pi * Dmid / MID.s2;
  nBlades2 = floor(nBlades2);
  if rem(nBlades2,2)
      nBlades2 = nBlades2 + 1;
  end

  MID.s2 = pi * Dmid / nBlades2;
  TIP.s2 = pi * Dtip / nBlades2;
  HUB.s2 = pi * Dhub / nBlades2;

  %MID.sigma1 IMPOSED;
  TIP.sigma2 =  TIP.c2 / TIP.s2;
  HUB.sigma2 = HUB.c2 / HUB.s2;

  %% Camber angle --> need iteration
  % epsilon = deltaBeta = theta + i - delta
  HUB.epsilon2 = HUB.deltaBeta2;
  MID.epsilon2 = MID.deltaBeta2;
  TIP.epsilon2 = TIP.deltaBeta2;

  %guess
  HUB.theta2 = HUB.epsilon2;
  MID.theta2 = MID.epsilon2;
  TIP.theta2 = TIP.epsilon2;

  nBlades=nBlades2;

  HUB.Re2 = Re;
  MID.Re2 = Re;
  TIP.Re2 = Re;

end

end
