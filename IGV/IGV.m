% IGV

% Off-Design Mass Flow Rate
mdot=90;

% New Axial Velocity from Section
HUBfp.v1a = mdot / (S * rho1);
MIDfp.v1a = mdot / (S * rho1);
TIPfp.v1a = mdot / (S * rho1);

HUBfp.w1t=tand(HUB.beta1)*MIDfp.v1a;
MIDfp.w1t=tand(MID.beta1)*MIDfp.v1a;
TIPfp.w1t=tand(TIP.beta1)*MIDfp.v1a;

HUBfp.v1t=HUBfp.w1t+U1(1);
MIDfp.v1t=MIDfp.w1t+U1(2);
TIPfp.v1t=TIPfp.w1t+U1(3);

HUBfp.v1=sqrt(HUBfp.v1t^2+MIDfp.v1a^2);
MIDfp.v1=sqrt(MIDfp.v1t^2+MIDfp.v1a^2);
TIPfp.v1=sqrt(TIPfp.v1t^2+MIDfp.v1a^2);

HUBfp.alfa1=atand(HUBfp.v1t/MIDfp.v1a);
alfa1=atand(MIDfp.v1t/MIDfp.v1a);
TIPfp.alfa1=atand(TIPfp.v1t/MIDfp.v1a);

% Keep fixed geometrical angles
HUB.beta2G = HUB.beta2 - HUB.delta1;
MID.beta2G = HUB.beta2 - HUB.delta1;
TIP.beta2G = HUB.beta2 - HUB.delta1;

HUB.beta3G = HUB.beta3 - HUB.i_opt2;
MID.beta3G = HUB.beta3 - HUB.i_opt2;
TIP.beta3G = HUB.beta3 - HUB.i_opt2;

HUB.beta4G = HUB.beta4 - HUB.delta2;
MID.beta4G = HUB.beta4 - HUB.delta2;
TIP.beta4G = HUB.beta4 - HUB.delta2;

% Estimate New Relative angles
HUBfp.beta2 = HUB.beta2;
MIDfp.beta2 = MID.beta2;
TIPfp.beta2 = TIP.beta2;

iter_fp_rot1 = 0;
maxiter = 100;
tol = 1e-3;
err_axvel_rot1 = 0;

hub_v2a_old= MIDfp.v1a;
mid_v2a_old= MIDfp.v1a;
tip_v2a_old= MIDfp.v1a;

while err_axvel_rot1 > tol && iter_fp_rot1 < maxiter || iter_fp_rot1 == 0;

  HUBfp.v2a = hub_v2a_old;
  MIDfp.v2a = mid_v2a_old;
  TIPfp.v2a = tip_v2a_old;

  HUBfp.w2t = tand(HUBfp.beta2) * HUBfp.v2a;
  MIDfp.w2t = tand(MIDfp.beta2) * MIDfp.v2a;
  TIPfp.w2t = tand(TIPfp.beta2) * TIPfp.v2a;

  HUBfp.v2t = U1(1) - HUBfp.w2t;
  MIDfp.v2t = U1(2) - MIDfp.w2t;
  TIPfp.v2t = U1(3) - TIPfp.w2t;

  HUBfp.leul1 = U1(1) * (HUBfp.v2t - HUBfp.v1t);
  MIDfp.leul1 = U1(2) * (MIDfp.v2t - MIDfp.v1t);
  TIPfp.leul1 = U1(3) * (TIP.v2t - TIPfp.v1t);

  HUBfp.v1 = sqrt(HUBfp.v1a^2 + HUBfp.v1t^2);
  MIDfp.v1 = sqrt(MIDfp.v1a^2 + MIDfp.v1t^2);
  TIPfp.v1 = sqrt(TIPfp.v1a^2 + TIPfp.v1t^2);

  HUBfp.v2 = sqrt(HUBfp.v2a^2 + HUBfp.v2t^2);
  MIDfp.v2 = sqrt(MIDfp.v2a^2 + MIDfp.v2t^2);
  TIPfp.v2 = sqrt(TIPfp.v2a^2 + TIPfp.v2t^2);

  HUBfp.Dhis1 = HUBfp.leul1 * etaTT - (HUBfp.v2^2 - HUBfp.v1^2) / 2;
  MIDfp.Dhis1 = MIDfp.leul1 * etaTT - (MIDfp.v2^2 - MIDfp.v1^2) / 2;
  TIPfp.Dhis1 = TIPfp.leul1 * etaTT - (TIPfp.v2^2 - TIPfp.v1^2) / 2;

  HUBfp.T2 = T1 + HUBfp.Dhis1 / Cp / etaTT;
  MIDfp.T2 = T1 + MIDfp.Dhis1 / Cp / etaTT;
  TIPfp.T2 = T1 + TIPfp.Dhis1 / Cp / etaTT;

  HUBfp.Beta1 = (1 + HUBfp.Dhis1 / (Cp * T1))^(gamma / (gamma - 1));
  MIDfp.Beta1 = (1 + MIDfp.Dhis1 / (Cp * T1))^(gamma / (gamma - 1));
  TIPfp.Beta1 = (1 + TIPfp.Dhis1 / (Cp * T1))^(gamma / (gamma - 1));

  HUBfp.P2 = HUBfp.Beta1 * P1;
  MIDfp.P2 = MIDfp.Beta1 * P1;
  TIPfp.P2 = TIPfp.Beta1 * P1;

  HUBfp.rho2 = HUBfp.P2 / (R * HUBfp.T2);
  MIDfp.rho2 = MIDfp.P2 / (R * MIDfp.T2);
  TIPfp.rho2 = TIPfp.P2 / (R * TIPfp.T2);


  %Update v2a
  HUBfp.v2a = mdot / (S * HUBfp.rho2);
  MIDfp.v2a = mdot / (S * MIDfp.rho2);
  TIPfp.v2a = mdot / (S * TIPfp.rho2);


  err_axvel_rot1(1) = abs( hub_v2a_old - HUBfp.v2a) / hub_v2a_old;
  err_axvel_rot1(2) = abs( mid_v2a_old - MIDfp.v2a) / mid_v2a_old;
  err_axvel_rot1(3) = abs( tip_v2a_old - TIPfp.v2a) / tip_v2a_old;

  err_axvel_rot1 = max(err_axvel_rot1);

  hub_v2a_old = HUBfp.v2a;
  mid_v2a_old = MIDfp.v2a;
  tip_v2a_old = TIPfp.v2a;

  iter_fp_rot1 = iter_fp_rot1 + 1;

end


% Estimate New Relative anglesatand(w3t/w3a)


HUBfp.w3a = HUBfp.v2a;
MIDfp.w3a = MIDfp.v2a;
TIPfp.w3a = TIPfp.v2a;

HUBfp.w3t = HUBfp.v2t - U2(1);
MIDfp.w3t = MIDfp.v2t - U2(2);
TIPfp.w3t = TIPfp.v2t - U2(3);

HUBfp.beta3 = atand(HUBfp.w3t / HUBfp.w3a);
MIDfp.beta3 = atand(MIDfp.w3t / MIDfp.w3a);
TIPfp.beta3 = atand(TIPfp.w3t / TIPfp.w3a);

HUBfp.i2 = HUBfp.beta3 - HUBfp.beta3G;
MIDfp.i2 = MIDfp.beta3 - MIDfp.beta3G;
TIPfp.i2 = TIPfp.beta3 - TIPfp.beta3G;

HUBfp.epsilon2 = HUBfp.i_opt2 + ( - HUBfp.beta4 + HUBfp.beta4G) + HUBfp.theta2;
MIDfp.epsilon2 = MIDfp.i_opt2 + ( - MIDfp.beta4 + MIDfp.beta4G) + MIDfp.theta2;
TIPfp.epsilon2 = TIPfp.i_opt2 + ( - TIPfp.beta4 + TIPfp.beta4G) + TIPfp.theta2;

HUBfp.beta4 = HUBfp.epsilon2 + HUBfp.beta3;
MIDfp.beta4 = MIDfp.epsilon2 + MIDfp.beta3;
TIPfp.beta4 = TIPfp.epsilon2 + TIPfp.beta3;


iter_fp_rot2 = 0;
maxiter = 100;
tol = 1e-3;
err_axvel_rot2 = 0;

hub_v4a_old= MIDfp.v2a;
mid_v4a_old= MIDfp.v2a;
tip_v4a_old= MIDfp.v2a;

while err_axvel_rot2 > tol && iter_fp_rot2 < maxiter || iter_fp_rot2 == 0;

  HUBfp.v4a = hub_v4a_old;
  MIDfp.v4a = mid_v4a_old;
  TIPfp.v4a = tip_v4a_old;

  HUBfp.w4t = tand(HUBfp.beta4) * HUBfp.v4a;
  MIDfp.w4t = tand(MIDfp.beta4) * MIDfp.v4a;
  TIPfp.w4t = tand(TIPfp.beta4) * TIPfp.v4a;

  HUBfp.v4t = U2(1) - HUBfp.w4t;
  MIDfp.v4t = U2(2) - MIDfp.w4t;
  TIPfp.v4t = U2(3) - TIPfp.w4t;

  HUBfp.leul2 = U2(1) * (HUBfp.v4t - HUBfp.v2t);
  MIDfp.leul2 = U2(2) * (MIDfp.v4t - MIDfp.v2t);
  TIPfp.leul2 = U2(3) * (TIP.v4t - TIPfp.v2t);

  HUBfp.v2 = sqrt(HUBfp.v2a^2 + HUBfp.v2t^2);
  MIDfp.v2 = sqrt(MIDfp.v2a^2 + MIDfp.v2t^2);
  TIPfp.v2 = sqrt(TIPfp.v2a^2 + TIPfp.v2t^2);

  HUBfp.v4 = sqrt(HUBfp.v4a^2 + HUBfp.v4t^2);
  MIDfp.v4 = sqrt(MIDfp.v4a^2 + MIDfp.v4t^2);
  TIPfp.v4 = sqrt(TIPfp.v4a^2 + TIPfp.v4t^2);

  HUBfp.Dhis2 = HUBfp.leul2 * etaTT - (HUBfp.v4^2 - HUBfp.v2^2) / 2;
  MIDfp.Dhis2 = MIDfp.leul2 * etaTT - (MIDfp.v4^2 - MIDfp.v2^2) / 2;
  TIPfp.Dhis2 = TIPfp.leul2 * etaTT - (TIPfp.v4^2 - TIPfp.v2^2) / 2;

  HUBfp.T4 = HUBfp.T2 + HUBfp.Dhis2 / Cp / etaTT;
  MIDfp.T4 = MIDfp.T2 + MIDfp.Dhis2 / Cp / etaTT;
  TIPfp.T4 = TIPfp.T2 + TIPfp.Dhis2 / Cp / etaTT;

  HUBfp.Beta2 = (1 + HUBfp.Dhis2 / (Cp * HUBfp.T2))^(gamma / (gamma - 1));
  MIDfp.Beta2 = (1 + MIDfp.Dhis2 / (Cp * MIDfp.T2))^(gamma / (gamma - 1));
  TIPfp.Beta2 = (1 + TIPfp.Dhis2 / (Cp * TIPfp.T2))^(gamma / (gamma - 1));

  HUBfp.P4 = HUBfp.Beta2 * HUBfp.P2;
  MIDfp.P4 = MIDfp.Beta2 * MIDfp.P2;
  TIPfp.P4 = TIPfp.Beta2 * TIPfp.P2;

  HUBfp.rho4 = HUBfp.P4 / (R * HUBfp.T4);
  MIDfp.rho4 = MIDfp.P4 / (R * MIDfp.T4);
  TIPfp.rho4 = TIPfp.P4 / (R * TIPfp.T4);


  %Update v2a
  HUBfp.v4a = mdot / (S * HUBfp.rho4);
  MIDfp.v4a = mdot / (S * MIDfp.rho4);
  TIPfp.v4a = mdot / (S * TIPfp.rho4);


  err_axvel_rot2(1) = abs( hub_v4a_old - HUBfp.v4a) / hub_v4a_old;
  err_axvel_rot2(2) = abs( mid_v4a_old - MIDfp.v4a) / mid_v4a_old;
  err_axvel_rot2(3) = abs( tip_v4a_old - TIPfp.v4a) / tip_v4a_old;

  err_axvel_rot2 = max(err_axvel_rot2);

  hub_v4a_old = HUBfp.v4a;
  mid_v4a_old = MIDfp.v4a;
  tip_v4a_old = TIPfp.v4a;

  iter_fp_rot2 = iter_fp_rot2 + 1;

end




HUBfp.v2t = U1(1) - HUBfp.w2t;
MIDfp.v2t = U1(2) - MIDfp.w2t;
TIPfp.v2t = U1(3) - TIPfp.w2t;

HUBfp.w3t=HUBfp.v2t-U2(1);
MIDfp.w3t=MIDfp.v2t-U2(2);
TIPfp.w3t=TIPfp.v2t-U2(3);



% % % % % % %Guess Beta_FP and try to match the reduced mass flow rate
% % % % % % Beta_new=1.4;
% % % % % % i_mdotfp = 0;
% % % % % % err_deltaBeta = 0;
% % % % % % maxiter = 100000;
% % % % % % tol = 1e-2;
% % % % % % MIDfp.v2a=MIDfp.v1a;
% % % % % % MIDfp.v4a=MIDfp.v1a;
% % % % % % mdot_iter=mdot;
% % % % % %
% % % % % % while err_deltaBeta > tol && i_mdotfp < maxiter || i_mdotfp == 0
% % % % % % BetaTot = Beta_new;
% % % % % %
% % % % % % Dhis=Cp*T1*(BetaTot^((gamma-1)/gamma)-1);
% % % % % %
% % % % % %
% % % % % % [MIDfp leulTotfp leul1fp leul2fp]=velocityTriangles(mdot_iter,alfa1,MIDfp.v1, MIDfp.v1a, MIDfp.v2a, MIDfp.v4a,S,P1,T1,U1(2),U2(2), etaTT,[0 1 0],Dhis,work1);
% % % % % %
% % % % % %
% % % % % % err_deltaBeta = abs(MID.deltaBeta1-MIDfp.deltaBeta1);
% % % % % %
% % % % % % i_mdotfp = i_mdotfp + 1;
% % % % % % Beta_new = BetaTot + 0.0001;
% % % % % % if Beta_new > 2
% % % % % %   disp('Non c''è soluzione')
% % % % % %   return
% % % % % %  end
% % % % % % end
