%IGV

% New Mass Flow Rate
mdot=90;

% Inlet Axial Velocity
HUBfp.v1a = mdot / (S * rho1);
MIDfp.v1a = mdot / (S * rho1);
TIPfp.v1a = mdot / (S * rho1);

HUBfp.v2a = HUBfp.v1a;
MIDfp.v2a = MIDfp.v1a;
TIPfp.v2a = TIPfp.v1a;

% Incidence (starting maximum value)
HUBfp.i1 = HUB.i_opt1;
MIDfp.i1 = MID.i_opt1;
TIPfp.i1 = TIP.i_opt1;

for iterIGV=1:100;

  HUBfp.i1 = HUB.i_opt1-0.1*iterIGV;
  MIDfp.i1 = MID.i_opt1-0.1*iterIGV;
  TIPfp.i1 = TIP.i_opt1-0.1*iterIGV;
% beta1 from incidence
[HUBfp.beta1 MIDfp.beta1 TIPfp.beta1] = beta1fromi(HUB,MID,TIP,1);

[HUBfp.delta1 MIDfp.delta1 TIPfp.delta1] = deltaIGV1(HUBfp,MIDfp,TIPfp,HUB,MID,TIP,'DCA',thick1);

HUBfp.deltaBeta1 = HUBfp.i1 - HUBfp.delta1 + HUB.theta1;
MIDfp.deltaBeta1 = MIDfp.i1 - MIDfp.delta1 + MID.theta1;
TIPfp.deltaBeta1 = TIPfp.i1 - TIPfp.delta1 + TIP.theta1;

HUBfp.beta2 = HUBfp.beta1 + HUBfp.deltaBeta1;
MIDfp.beta2 = MIDfp.beta1 + MIDfp.deltaBeta1;
TIPfp.beta2 = TIPfp.beta1 + TIPfp.deltaBeta1;

%Velocities 1st Rotor
HUBfp.w1t = HUBfp.v1a * tand(HUBfp.beta1);
MIDfp.w1t = MIDfp.v1a * tand(HUBfp.beta1);
TIPfp.w1t = TIPfp.v1a * tand(HUBfp.beta1);

HUBfp.v1t = HUBfp.w1t + U1(1);
MIDfp.v1t = MIDfp.w1t + U1(2);
TIPfp.v1t = TIPfp.w1t + U1(3);

HUBfp.v1 = sqrt(HUBfp.v1t^2+HUBfp.v1a^2);
MIDfp.v1 = sqrt(MIDfp.v1t^2+MIDfp.v1a^2);
TIPfp.v1 = sqrt(TIPfp.v1t^2+TIPfp.v1a^2);

HUBfp.alfa1 = asind(HUBfp.v1t/HUBfp.v1);
MIDfp.alfa1 = asind(MIDfp.v1t/MIDfp.v1);
TIPfp.alfa1 = asind(TIPfp.v1t/TIPfp.v1);

%iterate on axial velocity
tol_IGV = 1e-3;
maxiter_IGV = 100;
err_vax_IGV = 0;
i_ax_IGV1 = 0;

while err_vax_IGV > tol_IGV && i_ax_IGV1 < maxiter_IGV || i_ax_IGV1 == 0

HUBfp.w2t = HUBfp.v2a * tand(HUBfp.beta2);
MIDfp.w2t = MIDfp.v2a * tand(HUBfp.beta2);
TIPfp.w2t = TIPfp.v2a * tand(HUBfp.beta2);

HUBfp.v2t = HUBfp.w2t + U1(1);
MIDfp.v2t = MIDfp.w2t + U1(2);
TIPfp.v2t = TIPfp.w2t + U1(3);

HUBfp.v2 = sqrt(HUBfp.v2t^2+HUBfp.v2a^2);
MIDfp.v2 = sqrt(MIDfp.v2t^2+MIDfp.v2a^2);
TIPfp.v2 = sqrt(TIPfp.v2t^2+TIPfp.v2a^2);

HUBfp.alfa2 = asind(HUBfp.v2t/HUBfp.v2);
MIDfp.alfa2 = asind(MIDfp.v2t/MIDfp.v2);
TIPfp.alfa2 = asind(TIPfp.v2t/TIPfp.v2);

HUBfp.leul1 = U1(1)*(HUBfp.v2t-HUBfp.v1t);
MIDfp.leul1 = U1(2)*(MIDfp.v2t-MIDfp.v1t);
TIPfp.leul1 = U1(3)*(TIPfp.v2t-TIPfp.v1t);

HUBfp.Dhis1=HUBfp.leul1*eta1-(HUBfp.v2^2-HUBfp.v1^2)/2;
MIDfp.Dhis1=MIDfp.leul1*eta1-(MIDfp.v2^2-MIDfp.v1^2)/2;
TIPfp.Dhis1=TIPfp.leul1*eta1-(TIPfp.v2^2-TIPfp.v1^2)/2;

HUBfp.T2=T1+HUBfp.Dhis1/Cp/eta1;
MIDfp.T2=T1+MIDfp.Dhis1/Cp/eta1;
TIPfp.T2=T1+TIPfp.Dhis1/Cp/eta1;

HUBfp.Beta1=(1+HUBfp.Dhis1/(Cp*T1))^(gamma/(gamma-1));
MIDfp.Beta1=(1+MIDfp.Dhis1/(Cp*T1))^(gamma/(gamma-1));
TIPfp.Beta1=(1+TIPfp.Dhis1/(Cp*T1))^(gamma/(gamma-1));

HUBfp.P2=HUBfp.Beta1*P1;
MIDfp.P2=MIDfp.Beta1*P1;
TIPfp.P2=TIPfp.Beta1*P1;

HUBfp.rho2=HUBfp.P2/(R*HUBfp.T2);
MIDfp.rho2=MIDfp.P2/(R*MIDfp.T2);
TIPfp.rho2=TIPfp.P2/(R*TIPfp.T2);

v2a_hubnew=mdot/(S*HUBfp.rho2);
v2a_midnew=mdot/(S*MIDfp.rho2);
v2a_tipnew=mdot/(S*TIPfp.rho2);


err_vax_IGV(1) = abs(HUBfp.v2a - v2a_hubnew)/HUBfp.v2a;
err_vax_IGV(2) = abs(MIDfp.v2a - v2a_midnew)/MIDfp.v2a;
err_vax_IGV(3) = abs(TIPfp.v2a - v2a_tipnew)/TIPfp.v2a;

HUBfp.v2a = v2a_hubnew;
MIDfp.v2a = v2a_midnew;
TIPfp.v2a = v2a_tipnew;

err_vax_IGV = max(err_vax_IGV);

i_ax_IGV1 = i_ax_IGV1 + 1;

end

clearvars v2a_hubnew v2a_midnew v2a_tipnew


%velocities 2nd Rotor

HUBfp.w3t = HUBfp.v2t - U2(1);
MIDfp.w3t = MIDfp.v2t - U2(2);
TIPfp.w3t = TIPfp.v2t - U2(3);

HUBfp.beta3 = atand(HUBfp.w3t/HUBfp.v2a);
MIDfp.beta3 = atand(MIDfp.w3t/MIDfp.v2a);
TIPfp.beta3 = atand(TIPfp.w3t/TIPfp.v2a);

%lieblein_rot2_IGV;
% LIeblein_rot2_IGV
i010_tip = lieblein_incidence010(TIPfp.beta2,TIP.sigma2);%4;
i010_mid = lieblein_incidence010(MIDfp.beta2,MID.sigma2);%4;
i010_hub = lieblein_incidence010(HUBfp.beta2,HUB.sigma2);%4.2;

%correction coefficient for the shape
Ki_shape = lieblein_icorrectionShape(profile{4});

%correction for the thickness
Ki_thick_tip = lieblein_icorrectionThick(thick2); %0.9;
Ki_thick_mid = lieblein_icorrectionThick(thick2); %0.9;
Ki_thick_hub = lieblein_icorrectionThick(thick2); %0.9;

% n coefficient
% do graph --> [n_lieblein] = lieblein_ncoeff(beta1,sigma)
n_lieblein_tip = lieblein_ncoeff(TIPfp.beta2,TIP.sigma2); %-0.26;
n_lieblein_mid = lieblein_ncoeff(MIDfp.beta2,MID.sigma2); %-0.21;
n_lieblein_hub = lieblein_ncoeff(HUBfp.beta2,HUB.sigma2); %-0.19;

i0_tip = Ki_thick_tip * Ki_shape * i010_tip;
i0_mid = Ki_thick_mid * Ki_shape * i010_mid;
i0_hub = Ki_thick_hub * Ki_shape * i010_hub;

HUBfp.i2 = i0_hub + n_lieblein_hub * HUB.theta2;
MIDfp.i2 = i0_mid + n_lieblein_mid * MID.theta2;
TIPfp.i2 = i0_tip + n_lieblein_tip * TIP.theta2;


[HUBfp.delta2 MIDfp.delta2 TIPfp.delta2] = deltaIGV2(HUBfp,MIDfp,TIPfp,HUB,MID,TIP,'DCA',thick2);

HUBfp.deltaBeta2 = HUBfp.i2 - HUBfp.delta2 + HUB.theta2;
MIDfp.deltaBeta2 = MIDfp.i2 - MIDfp.delta2 + MID.theta2;
TIPfp.deltaBeta2 = TIPfp.i2 - TIPfp.delta2 + TIP.theta2;

HUBfp.beta4 = HUBfp.beta3 + HUBfp.deltaBeta2;
MIDfp.beta4 = MIDfp.beta3 + MIDfp.deltaBeta2;
TIPfp.beta4 = TIPfp.beta3 + TIPfp.deltaBeta2;

% Iterate on axial velocity
HUBfp.v4a = HUBfp.v2a;
MIDfp.v4a = MIDfp.v2a;
TIPfp.v4a = TIPfp.v2a;


i_ax_IGV2 = 0;

while err_vax_IGV > tol_IGV && i_ax_IGV2 < maxiter_IGV || i_ax_IGV2 == 0

HUBfp.w4t = HUBfp.v4a * tand(HUBfp.beta4);
MIDfp.w4t = MIDfp.v4a * tand(HUBfp.beta4);
TIPfp.w4t = TIPfp.v4a * tand(HUBfp.beta4);

HUBfp.v4t = HUBfp.w4t + U2(1);
MIDfp.v4t = MIDfp.w4t + U2(2);
TIPfp.v4t = TIPfp.w4t + U2(3);

HUBfp.v4 = sqrt(HUBfp.v4t^2+HUBfp.v4a^2);
MIDfp.v4 = sqrt(MIDfp.v4t^2+MIDfp.v4a^2);
TIPfp.v4 = sqrt(TIPfp.v4t^2+TIPfp.v4a^2);

HUBfp.alfa4 = asind(HUBfp.v4t/HUBfp.v4);
MIDfp.alfa4 = asind(MIDfp.v4t/MIDfp.v4);
TIPfp.alfa4 = asind(TIPfp.v4t/TIPfp.v4);

HUBfp.leul2 = U2(1)*(HUBfp.v4t-HUBfp.v2t);
MIDfp.leul2 = U2(2)*(MIDfp.v4t-MIDfp.v2t);
TIPfp.leul2 = U2(3)*(TIPfp.v4t-TIPfp.v2t);

HUBfp.Dhis2=HUBfp.leul2*eta2-(HUBfp.v4^2-HUBfp.v2^2)/2;
MIDfp.Dhis2=MIDfp.leul2*eta2-(MIDfp.v4^2-MIDfp.v2^2)/2;
TIPfp.Dhis2=TIPfp.leul2*eta2-(TIPfp.v4^2-TIPfp.v2^2)/2;

HUBfp.T4=HUBfp.T2+HUBfp.Dhis2/Cp/eta2;
MIDfp.T4=MIDfp.T2+MIDfp.Dhis2/Cp/eta2;
TIPfp.T4=TIPfp.T2+TIPfp.Dhis2/Cp/eta2;

HUBfp.Beta2=(1+HUBfp.Dhis2/(Cp*HUBfp.T2))^(gamma/(gamma-1));
MIDfp.Beta2=(1+MIDfp.Dhis2/(Cp*MIDfp.T2))^(gamma/(gamma-1));
TIPfp.Beta2=(1+TIPfp.Dhis2/(Cp*TIPfp.T2))^(gamma/(gamma-1));

HUBfp.P4=HUBfp.Beta2*HUBfp.P2;
MIDfp.P4=MIDfp.Beta2*MIDfp.P2;
TIPfp.P4=TIPfp.Beta2*TIPfp.P2;

HUBfp.rho4=HUBfp.P4/(R*HUBfp.T4);
MIDfp.rho4=MIDfp.P4/(R*MIDfp.T4);
TIPfp.rho4=TIPfp.P4/(R*TIPfp.T4);

v4a_hubnew=mdot/(S*HUBfp.rho4);
v4a_midnew=mdot/(S*MIDfp.rho4);
v4a_tipnew=mdot/(S*TIPfp.rho4);

err_vax_IGV(1) = abs(HUBfp.v4a - v4a_hubnew)/HUBfp.v4a;
err_vax_IGV(2) = abs(MIDfp.v4a - v4a_midnew)/MIDfp.v4a;
err_vax_IGV(3) = abs(TIPfp.v4a - v4a_tipnew)/TIPfp.v4a;

HUBfp.v4a = v4a_hubnew;
MIDfp.v4a = v4a_midnew;
TIPfp.v4a = v4a_tipnew;

err_vax_IGV = max(err_vax_IGV);

i_ax_IGV2 = i_ax_IGV2 + 1;

end

clearvars v4a_hubnew v4a_midnew v4a_tipnew


%% BLADES, loading (choose solidity) and camber angle

%Howell correlation per calcolare il
iterHowell_IGV = 0;
maxiter_IGV = 100;
tol = 1e-3;

sigmaLimit = 2.5;
sigma=0.625:0.01:sigmaLimit;


%while changeSolidity < 0.01 && changeSolidity > -0.01 || sigma > sigmaLimit && iterHowell < maxiter || iterHowell == 0

for i = 1 : length(sigma)
[deltaBetaOpt1(1)] = howellCorrelation(HUBfp.alfa1,5e5,sigma(i));
[deltaBetaOpt1(2)] = howellCorrelation(MIDfp.alfa1,5e5,sigma(i));
[deltaBetaOpt1(3)] = howellCorrelation(TIPfp.alfa1,5e5,sigma(i));

changeSolidity1(1)=deltaBetaOpt1(1)-HUBfp.deltaBeta1;  %se <1 aumentare solidity
changeSolidity1(2)=deltaBetaOpt1(2)-MIDfp.deltaBeta1;  %se >1 diminuire solidity
changeSolidity1(3)=deltaBetaOpt1(3)-TIPfp.deltaBeta1;
changeSolidity(:,i) = [changeSolidity1'; changeSolidity2'];

end

[changeBest,index] = min(abs(changeSolidity)');
maxDeflAllowed = 3;
if changeBest(1) > maxDeflAllowed || ...
        changeBest(2) > maxDeflAllowed || ...
        changeBest(3) > maxDeflAllowed
    disp('howell_1_IGV non è ok')

    solidityBest = sigma(index);
    %return
else
solidityBest = sigma(index);

end

%% New Part
HUBfp.deltaAlfa = HUBfp.alfa1 - 0;
MIDfp.deltaAlfa = MIDfp.alfa1 - 0;
TIPfp.deltaAlfa = TIPfp.alfa1 - 0;


%number of blades
[HUBfp, MIDfp, TIPfp, nBladesIGV] = nBladesFUNC_IGV(HUBfp,MIDfp,TIPfp,rho1,mu,Re1,solidityBest,Dhub,Dmid,Dtip);

%NEW Howell

[deltaBetaOpt1(1)] = howellCorrelation(HUB.beta2,5e5,HUB.sigma1);
[deltaBetaOpt1(2)] = howellCorrelation(MID.beta2,5e5,MID.sigma1);
[deltaBetaOpt1(3)] = howellCorrelation(TIP.beta2,5e5,TIP.sigma1);

changeSolidity1(1)=deltaBetaOpt1(1)-HUBfp.deltaBeta1;  %se <1 aumentare solidity
changeSolidity1(2)=deltaBetaOpt1(2)-MIDfp.deltaBeta1;  %se >1 diminuire solidity
changeSolidity1(3)=deltaBetaOpt1(3)-TIPfp.deltaBeta1;

for i = 1 : length(sigma)
changeSolidity(:,i) = [changeSolidity1'; changeSolidity2'];
end

[changeBest,~] = min(abs(changeSolidity)')
maxDeflAllowed = 3;
if changeBest(1) > maxDeflAllowed || ...
        changeBest(2) > maxDeflAllowed || ...
        changeBest(3) > maxDeflAllowed
    disp('howell_2_IGV non è ok')
end


% CHOICE!!!!!!!!!!!!!!!!
%thick1 = 0.1;      %thickness 10% of the chord
%thick2 = 0.1;

%Iteration on CAMBER ANGLE
i = 0;
maxiter = 100;
tol = 1e-3;
err = 0;

while err > tol && i < maxiter || i == 0

%% Camber angle --> need iteration
% epsilon = deltaBeta = theta + i - delta

%START WITH NACA PROFILES, link between CL and the angle
% FOR NACA 65
%Cl = @(theta) theta / 25;
% FOR DCA
Cl = @(theta) tand(theta/4) / 0.1103;

%theta is the camber angle
HUBfp.Cl_IGV = Cl(HUBfp.theta_IGV);
MIDfp.Cl_IGV = Cl(MIDfp.theta_IGV);
TIPfp.Cl_IGV = Cl(TIPfp.theta_IGV);

%Re = @(rho,u,d) rho * u * d / mu;
%
%Re_hub = Re((rho1+rho2_hub)/2,(w1_hub + w2_hub)/2,HUB.c1);
%Re_mid = Re((rho1+rho2_mid)/2,(w1_mid + w2_mid)/2,c);
%Re_tip = Re((rho1+rho2_tip)/2,(w1_tip + w2_tip)/2,c);

lieblein_IGV;

end

%keyboard
checkLoading_IGV;

%stagger angle
HUBfp.gamma_IGV = 0 + HUBfp.i_IGV + HUBfp.theta_IGV / 2;
MIDfp.gamma_IGV = 0 + MIDfp.i_IGV + MIDfp.theta_IGV / 2;
TIPfp.gamma_IGV = 0 + TIPfp.i_IGV + TIPfp.theta_IGV / 2;

%% Losses
%% Profile Losses
%profileLosses;
saveGeometryfp

[HUBfp.Zprofile1, MIDfp.Zprofile1, TIPfp.Zprofile1] = profileLosses_rot1('traupel',profile, HUBfp, MIDfp, TIPfp);
[HUBfp.Zprofile2, MIDfp.Zprofile2, TIPfp.Zprofile2] = profileLosses_rot2('traupel',profile, HUBfp, MIDfp, TIPfp);

%% Tip clearance / Leakage losses
deltaC = 2e-3; % tip leakage gap
[deltaEta_leakage1_fp] = tipClearanceLoss(mdot,deltaC,b,TIPfp.c1,U1(2),MIDfp.v1a,MIDfp.leul1,TIPfp.beta1,TIPfp.beta2,eta1,S,rho1);
[deltaEta_leakage2_fp] = tipClearanceLoss(mdot,deltaC,b,TIPfp.c2,U2(2),MIDfp.v2a,MIDfp.leul2,TIPfp.beta3,TIPfp.beta4,eta2,S,MIDfp.rho2);

%% Annulus Losses
%% Secondary Losses / Endwall

[zetaEW1_fp] = endwallLosses(HUBfp.Zprofile1,HUBfp.P2,HUBfp.Pt2,HUBfp.M2,HUBfp.deltaBeta1,rho1,HUBfp.rho2,HUBfp.v1,HUBfp.v2,HUBfp.alfa2,t_TE1,b,S);
[zetaEW2_fp] = endwallLosses(HUBfp.Zprofile2,HUBfp.P4,HUBfp.Pt4,HUBfp.M4,HUBfp.deltaBeta2,HUBfp.rho2,HUBfp.rho4,HUBfp.v2,HUBfp.v4,HUBfp.alfa4,t_TE2,b,S);

%% Overall Losses Combination (this fromula from Y to Zeta then to deltaEta)
% Yprofile --> ZetaP --> deltaEta_profile
[HUBfp.deltaEta_profile1] = deltaEta_calc('Z',[],HUBfp.Zprofile1,HUBfp.Mw2,P1,HUBfp.P2,T1,HUBfp.Beta1,Pt1,v1);
[MIDfp.deltaEta_profile1] = deltaEta_calc('Z',[],MIDfp.Zprofile1,MIDfp.Mw2,P1,MIDfp.P2,T1,MIDfp.Beta1,Pt1,v1);
[TIPfp.deltaEta_profile1] = deltaEta_calc('Z',[],TIPfp.Zprofile1,TIPfp.Mw2,P1,TIPfp.P2,T1,TIPfp.Beta1,Pt1,v1);

[HUBfp.deltaEta_profile2] = deltaEta_calc('Z',[],HUBfp.Zprofile2,HUBfp.Mw4,HUBfp.P2,HUBfp.P4,HUBfp.T2,HUBfp.Beta2,HUBfp.Pt2,HUBfp.v2);
[MIDfp.deltaEta_profile2] = deltaEta_calc('Z',[],MIDfp.Zprofile2,MIDfp.Mw4,MIDfp.P2,MIDfp.P4,MIDfp.T2,MIDfp.Beta2,MIDfp.Pt2,MIDfp.v2);
[TIPfp.deltaEta_profile2] = deltaEta_calc('Z',[],TIPfp.Zprofile2,TIPfp.Mw4,TIPfp.P2,TIPfp.P4,TIPfp.T2,TIPfp.Beta2,TIPfp.Pt2,TIPfp.v2);

% zetaEW --> deltaEta_endwall
[deltaEta_endwall1_fp] = deltaEta_calc('Z',[],zetaEW1_fp,HUBfp.M2,P1,HUBfp.P2,T1,HUBfp.Beta1,Pt1,HUBfp.v1);
[deltaEta_endwall2_fp] = deltaEta_calc('Z',[],zetaEW2_fp,HUBfp.M4,HUBfp.P2,HUBfp.P4,HUBfp.T2,HUBfp.Beta2,HUBfp.Pt2,HUBfp.v2);

eta1fp = 1-((HUBfp.deltaEta_profile1 + MIDfp.deltaEta_profile1 + TIPfp.deltaEta_profile1)/3 + deltaEta_endwall1_fp + deltaEta_leakage1_fp);
eta2fp = 1-((HUBfp.deltaEta_profile2 + MIDfp.deltaEta_profile2 + TIPfp.deltaEta_profile2)/3 + deltaEta_endwall2_fp + deltaEta_leakage2_fp);
etaTOTfp(iterIGV) = (( (MIDfp.Beta1*MIDfp.Beta2)^(GAMMA) - 1)*T1*eta1fp*eta2fp ) / (   T1*eta2fp*(MIDfp.Beta1^GAMMA-1) +  MIDfp.T2*eta1fp*(MIDfp.Beta2^GAMMA-1)   ) - deltaEta_diskFriction;

[deltaBetaOpt1(1)] = howellCorrelation(HUBfp.beta2,5e5,HUBfp.sigma1);
[deltaBetaOpt1(2)] = howellCorrelation(MIDfp.beta2,5e5,MIDfp.sigma1);
[deltaBetaOpt1(3)] = howellCorrelation(TIPfp.beta2,5e5,TIPfp.sigma1);
[deltaBetaOpt2(1)] = howellCorrelation(HUBfp.beta4,5e5,HUBfp.sigma2);
[deltaBetaOpt2(2)] = howellCorrelation(MIDfp.beta4,5e5,MIDfp.sigma2);
[deltaBetaOpt2(3)] = howellCorrelation(TIPfp.beta4,5e5,TIPfp.sigma2);

changeSolidity1(1)=deltaBetaOpt1(1)-HUB.deltaBeta1;  %se <1 aumentare solidity
changeSolidity1(2)=deltaBetaOpt1(2)-MID.deltaBeta1;  %se >1 diminuire solidity
changeSolidity1(3)=deltaBetaOpt1(3)-TIP.deltaBeta1;
changeSolidity2(1)=deltaBetaOpt2(1)+HUB.deltaBeta2;
changeSolidity2(2)=deltaBetaOpt2(2)+MID.deltaBeta2;
changeSolidity2(3)=deltaBetaOpt2(3)+TIP.deltaBeta2;

for i = 1 : length(sigma)
changeSolidity(:,i) = [changeSolidity1'; changeSolidity2'];
end

[changeBest,~] = min(abs(changeSolidity)');
if changeBest(1) > 4 || ...
        changeBest(2) > maxDeflAllowed || ...
        changeBest(3) > maxDeflAllowed || ...
        changeBest(4) > maxDeflAllowed || ...
        changeBest(5) > maxDeflAllowed || ...
        changeBest(6) > maxDeflAllowed
    etaTOTfp(iterIGV)=0.1;
end

end

[best_etaTOT iterBestIGV]= max(etaTOTfp);
