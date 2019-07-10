function [hub_Yprofile1 mid_Yprofile1 tip_Yprofile1] = profileLosses_rot1(name,profile, HUB, MID, TIP)

% profile needs to be a string containing the 3 profile names.
% According to Lieblein Correlation
hub_Yprofile1 = profileLossesLieblein(HUB.theta_c1,HUB.sigma1,HUB.beta1,HUB.beta2);
mid_Yprofile1 = profileLossesLieblein(MID.theta_c1,MID.sigma1,MID.beta1,MID.beta2);
tip_Yprofile1 = profileLossesLieblein(TIP.theta_c1,TIP.sigma1,TIP.beta1,TIP.beta2);

if strcmp(name,'leiblein')
%% Reynolds Correction
if Re > 3e5
  correction_Re = 0;
else
  %Grafico, non c'è quindi fisso
  correction_Re = 0.05;
end

%% Mach Correction
if Ma < 0.7
  correction_Ma = 0.015;
else
  %Grafico, non c'è ancora quindi fisso
end

%% Thickness Correction for Re

%% Thickness Correction for Ma

Yprofile = Yprofile - correction;

elseif strcmp(name,'traupel')
%% TRAUPEL
Zeta0_hub = hub_Yprofile1;
Zeta0_mid = mid_Yprofile1;
Zeta0_tip = tip_Yprofile1;


%% Re Correction (Traupel did calculation for Re = 2e5)
if HUB.Re1 > 4e4 && MID.Re1 > 4e4 && TIP.Re1 > 4e4

  CHI_Re_hub = (2e5 / HUB.Re1)^0.2;
  CHI_Re_mid = (2e5 / MID.Re1)^0.2;
  CHI_Re_tip = (2e5 / TIP.Re1)^0.2;

else

  CHI_Re_hub = 1.38 * (4e4 / HUB.Re1)^0.5;
  CHI_Re_mid = 1.38 * (4e4 / MID.Re1)^0.5;
  CHI_Re_tip = 1.38 * (4e4 / TIP.Re1)^0.5;

end

%% Ma Correction
flag = 1;

for i = 1 : 3

if strcmp(profile{i},'naca65')

  if i == 1
  [CHI_Ma_hub] = traupel_Ma_naca65(HUB.Mw1,HUB.sigma1,flag);
  elseif i == 2
  [CHI_Ma_mid] = traupel_Ma_naca65(MID.Mw1,MID.sigma1,flag);
  elseif i == 3
  [CHI_Ma_tip] = traupel_Ma_naca65(TIP.Mw1,TIP.sigma1,flag);
  end

elseif strcmp(profile{i},'DCA')

  if i == 1
    [CHI_Ma_hub] = traupel_Ma_transonic(HUB.Mw1,flag);
  elseif i == 2
    [CHI_Ma_mid] = traupel_Ma_transonic(MID.Mw1,flag);
  elseif i == 3
    [CHI_Ma_tip] = traupel_Ma_transonic(TIP.Mw1,flag);
  end

end

end

hub_Yprofile1 = CHI_Re_hub * CHI_Ma_hub * Zeta0_hub;
mid_Yprofile1 = CHI_Re_mid * CHI_Ma_mid * Zeta0_mid;
tip_Yprofile1 = CHI_Re_tip * CHI_Ma_tip * Zeta0_tip;

% Scaling to obtain absolute values from percentage
hub_Yprofile1 = hub_Yprofile1 * 0.01;
mid_Yprofile1 = mid_Yprofile1 * 0.01;
tip_Yprofile1 = tip_Yprofile1 * 0.01;

else
  disp('No profile losses criterion selected')
  return
end


end
