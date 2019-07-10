function [HUB, MID, TIP, nBlades] = nBladesFUNC_keepSIGMA(HUB,MID,TIP,rotor,rho,mu,Re,solidityBest,Dhub,Dmid,Dtip)

assignSolidityStruct;


if rotor==1

  nBlades1 = pi * Dmid / MID.sigma1 * rho * MID.w1 / (Re * mu);
  nBlades1 = floor(nBlades1);
  if rem(nBlades1,2)
      nBlades1 = nBlades1 + 1;
  end

  HUB.c1 = pi * Dhub / (nBlades1 * HUB.sigma1);
  MID.c1 = pi * Dmid / (nBlades1 * MID.sigma1);
  TIP.c1 = pi * Dtip / (nBlades1 * TIP.sigma1);

  HUB.s1 = HUB.c1 / HUB.sigma1;
  MID.s1 = MID.c1 / MID.sigma1;
  TIP.s1 = TIP.c1 / TIP.sigma1;

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

  %Check Re
  HUB.Re1 = rho * HUB.w1 * HUB.c1 / mu;
  MID.Re1 = Re;
  TIP.Re1 = rho * TIP.w1 * TIP.c1 / mu;

  %Check Re > 5e5
  if [HUB.Re1 MID.Re1 TIP.Re1] > 5e5
    disp('Re1 > 5e5 --> OK')
  end

elseif rotor==2
  nBlades2 = pi * Dmid / MID.sigma2 * MID.rho2 * MID.w2 / (Re * mu);
  nBlades2 = floor(nBlades2);
  if rem(nBlades2,2)
      nBlades2 = nBlades2 + 1;
  end

  HUB.c2 = pi * Dhub / (nBlades2 * HUB.sigma2);
  MID.c2 = pi * Dmid / (nBlades2 * MID.sigma2);
  TIP.c2 = pi * Dtip / (nBlades2 * TIP.sigma2);

  HUB.s2 = HUB.c2 / HUB.sigma2;
  MID.s2 = MID.c2 / MID.sigma2;
  TIP.s2 = TIP.c2 / TIP.sigma2;

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

  %Check Re
  HUB.Re2 = HUB.rho2 * HUB.w3 * HUB.c2 / mu;
  MID.Re2 = Re;
  TIP.Re2 = TIP.rho2 * TIP.w3 * TIP.c2 / mu;

  %Check Re > 5e5
  if [HUB.Re2 MID.Re2 TIP.Re2] > 5e5
    disp('Re2 > 5e5 --> OK')
  end

end



end
