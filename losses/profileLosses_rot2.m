function [hub_Zprofile2 mid_Zprofile2 tip_Zprofile2] = profileLosses_rot2(name,profile, HUB, MID, TIP)

% profile needs to be a string containing the 3 profile names.
% According to Lieblein Correlation
hub_Yprofile2 = profileLossesLieblein(HUB.theta_c2,HUB.sigma2,HUB.beta2,HUB.beta4);
mid_Yprofile2 = profileLossesLieblein(MID.theta_c2,MID.sigma2,MID.beta2,MID.beta4);
tip_Yprofile2 = profileLossesLieblein(TIP.theta_c2,TIP.sigma2,TIP.beta2,TIP.beta4);

%elseif strcmp(name,'traupel')
%% TRAUPEL
Zeta0_hub = hub_Yprofile2;
Zeta0_mid = mid_Yprofile2;
Zeta0_tip = tip_Yprofile2;


%% Re Correction (Traupel did calculation for Re = 2e5)
if HUB.Re2 > 4e4 && MID.Re2 > 4e4 && TIP.Re2 > 4e4

  CHI_Re_hub = (2e5 / HUB.Re2)^0.2;
  CHI_Re_mid = (2e5 / MID.Re2)^0.2;
  CHI_Re_tip = (2e5 / TIP.Re2)^0.2;

else

  CHI_Re_hub = 1.38 * (4e4 / HUB.Re2)^0.5;
  CHI_Re_mid = 1.38 * (4e4 / MID.Re2)^0.5;
  CHI_Re_tip = 1.38 * (4e4 / TIP.Re2)^0.5;

end

%% Ma Correction
flag = 1;

for i = 4 : 6

if strcmp(profile{i},'naca65')

  if i - 3 == 1
  [CHI_Ma_hub] = traupel_Ma_naca65(HUB.Mw3,HUB.sigma2,flag);
  elseif i - 3 == 2
  [CHI_Ma_mid] = traupel_Ma_naca65(MID.Mw3,MID.sigma2,flag);
  elseif i -3 == 3
  [CHI_Ma_tip] = traupel_Ma_naca65(TIP.Mw3,TIP.sigma2,flag);
  end

elseif strcmp(profile{i},'DCA')

  if i - 3 == 1
    [CHI_Ma_hub] = traupel_Ma_transonic(HUB.Mw3,flag);
  elseif i - 3 == 2
    [CHI_Ma_mid] = traupel_Ma_transonic(MID.Mw3,flag);
  elseif i - 3== 3
    [CHI_Ma_tip] = traupel_Ma_transonic(TIP.Mw3,flag);
  end

end

end

hub_Zprofile2 = CHI_Re_hub * CHI_Ma_hub * Zeta0_hub;
mid_Zprofile2 = CHI_Re_mid * CHI_Ma_mid * Zeta0_mid;
tip_Zprofile2 = CHI_Re_tip * CHI_Ma_tip * Zeta0_tip;


% Scaling to obtain absolute values from percentage
hub_Zprofile2 = hub_Zprofile2 * 0.01;
mid_Zprofile2 = mid_Zprofile2 * 0.01;
tip_Zprofile2 = tip_Zprofile2 * 0.01;


end
