clear all;
close all;
clc

%PATH
addPaths;

% DATA
mdot = 100;                 % Mass flow rate            [kg/s]
Pt1 = 100000;               % Inlet Total Pressure      [Pa]
Tt1 = 300;                  % Inlet Total Temperature   [K]
gamma = 1.4;                % Specific heat ratio
GAMMA = (gamma-1)/gamma;
R=287;
mu = 1.81e-5;
Cp = 1004.69;
BetaTot = 1.45;               % Required Compression Ratio

%Inlet guide vane is REQUIRED
%Free to choose rpm
%DESIGN WITH HIGHEST ACHIEVABLE EFFICIENCY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARAMETERS
Dtip=1;
eta1=0.76;
eta2=0.87;
load OPT5
opt=optimalSolution;

work1=opt(1);
rotRatio=opt(2);
n=opt(3);
Mw1_tip= opt(4);
Re1 = opt(5);
Re2 = opt(6);
alfa_mid=opt(7);
alfa_hub=opt(8);
alfa_tip=opt(9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega=n/60*2*pi;
U1_tip=omega*Dtip/2;

M1=2;
residual=100;

while residual>1e-3 && M1>0
    M1=M1-0.001;
    T1=Tt1/(1+(gamma-1)/2*M1^2);
    v1=sqrt(gamma*R*T1)*M1;
    w1t_tip=v1*sind(alfa_tip)-U1_tip;
    v1a_tip=v1*cosd(alfa_tip);
    w1_tip=sqrt(v1a_tip^2+w1t_tip^2);
    Mrel_tip_check=w1_tip/sqrt(gamma*R*T1);
    residual=abs(Mw1_tip-Mrel_tip_check);

end


v1a_mid=v1*cosd(alfa_mid);
v1a_hub=v1*cosd(alfa_hub);
clearvars U1_tip residual w1_tip w1t_tip

P1=Pt1*(T1/Tt1)^(gamma/(gamma-1));
rho1=P1/(R*T1);
S=mdot/(v1a_mid*rho1);
Dhub=@(Dhub) S-pi/4*(Dtip^2-Dhub^2);
Dhub=fzero(Dhub,Dtip);
b=(Dtip-Dhub)/2;
Dmid=(Dtip+Dhub)/2;

U1=[omega*Dhub/2 omega*Dmid/2 omega*Dtip/2];
U2=[-U1(1)*rotRatio -U1(2)*rotRatio -U1(3)*rotRatio]; %solo se b1=b2

Dhis=Cp*T1*(BetaTot^GAMMA-1);
v2a=v1a_mid; v4a=v1a_mid;


% Outer Iteration on Efficiency

i_eta = 0;
max_eta_iter = 100;
tol = 1e-2;
err_efficiency = 0;

while err_efficiency > tol && i_eta < max_eta_iter || i_eta == 0

[MID leulTot leul1 leul2]=velocityTriangles(mdot,alfa_mid,v1,v1a_mid,v2a,v4a,S,P1,T1,U1(2),U2(2),eta1,eta2,[0 1 0],Dhis,work1);
[HUB]=velocityTriangles(mdot,alfa_hub,v1,v1a_hub,MID.v2a,MID.v4a,S,P1,T1,U1(1),U2(1),eta1,eta2,[1 0 0],leulTot,work1,leul1,leul2);
[TIP]=velocityTriangles(mdot,alfa_tip,v1,v1a_tip,MID.v2a,MID.v4a,S,P1,T1,U1(3),U2(3),eta1,eta2,[0 0 1],leulTot,work1,leul1,leul2);

if HUB.Beta1 < 1
    disp('Beta < 1 @ hub, ERROR')
end

% Mach Number Check
Mrel_tip1=TIP.w1/sqrt(gamma*R*T1);
Mrel_tip2=TIP.w3/sqrt(gamma*R*TIP.T2);


%Howell correlation - Loading Criterion

%Provide a range of values. Find best solidity
sigmaLimit = 2.5;
sigma=0.625:0.01:sigmaLimit;

for i = 1 : length(sigma)
[deltaBetaOpt1(1)] = howellCorrelation(HUB.beta2,5e5,sigma(i));
[deltaBetaOpt1(2)] = howellCorrelation(MID.beta2,5e5,sigma(i));
[deltaBetaOpt1(3)] = howellCorrelation(TIP.beta2,5e5,sigma(i));
[deltaBetaOpt2(1)] = howellCorrelation(HUB.beta4,5e5,sigma(i));
[deltaBetaOpt2(2)] = howellCorrelation(MID.beta4,5e5,sigma(i));
[deltaBetaOpt2(3)] = howellCorrelation(TIP.beta4,5e5,sigma(i));

changeSolidity1(1)=deltaBetaOpt1(1)-HUB.deltaBeta1;  %if <1 increase solidity
changeSolidity1(2)=deltaBetaOpt1(2)-MID.deltaBeta1;  %if >1 decrease solidity
changeSolidity1(3)=deltaBetaOpt1(3)-TIP.deltaBeta1;
changeSolidity2(1)=deltaBetaOpt2(1)+HUB.deltaBeta2;
changeSolidity2(2)=deltaBetaOpt2(2)+MID.deltaBeta2;
changeSolidity2(3)=deltaBetaOpt2(3)+TIP.deltaBeta2;
changeSolidity(:,i) = [changeSolidity1'; changeSolidity2'];

end

%Find the position of the best solidities in the vector
[changeBest,index] = min(abs(changeSolidity)');

%Max difference allowed from Howell optimal deflection and actual deflection
maxDeflAllowed = 3;

if changeBest(1) > maxDeflAllowed || ...
        changeBest(2) > maxDeflAllowed || ...
        changeBest(3) > maxDeflAllowed || ...
        changeBest(4) > maxDeflAllowed || ...
        changeBest(5) > maxDeflAllowed || ...
        changeBest(6) > maxDeflAllowed
    disp('howell_1 non è ok')

    solidityBest = sigma(index);

else
    
% Assign Best solidity    
solidityBest = sigma(index);

end


% Number of blades
[HUB, MID, TIP, nBlades1] = nBladesFUNC(HUB,MID,TIP,1,rho1,mu,Re1,solidityBest,Dhub,Dmid,Dtip);
[HUB, MID, TIP, nBlades2] = nBladesFUNC(HUB,MID,TIP,2,[],mu,Re2,solidityBest,Dhub,Dmid,Dtip);

% Again Howell to check if the new solidity (has changed assigning the
% blades) is ok.

[deltaBetaOpt1(1)] = howellCorrelation(HUB.beta2,5e5,HUB.sigma1);
[deltaBetaOpt1(2)] = howellCorrelation(MID.beta2,5e5,MID.sigma1);
[deltaBetaOpt1(3)] = howellCorrelation(TIP.beta2,5e5,TIP.sigma1);
[deltaBetaOpt2(1)] = howellCorrelation(HUB.beta4,5e5,HUB.sigma2);
[deltaBetaOpt2(2)] = howellCorrelation(MID.beta4,5e5,MID.sigma2);
[deltaBetaOpt2(3)] = howellCorrelation(TIP.beta4,5e5,TIP.sigma2);

changeSolidity1(1)=deltaBetaOpt1(1)-HUB.deltaBeta1;
changeSolidity1(2)=deltaBetaOpt1(2)-MID.deltaBeta1;  
changeSolidity1(3)=deltaBetaOpt1(3)-TIP.deltaBeta1;
changeSolidity2(1)=deltaBetaOpt2(1)+HUB.deltaBeta2;
changeSolidity2(2)=deltaBetaOpt2(2)+MID.deltaBeta2;
changeSolidity2(3)=deltaBetaOpt2(3)+TIP.deltaBeta2;

for i = 1 : length(sigma)
changeSolidity(:,i) = [changeSolidity1'; changeSolidity2'];
end

[changeBest,~] = min(abs(changeSolidity)')

if changeBest(1) > maxDeflAllowed || ...
        changeBest(2) > maxDeflAllowed || ...
        changeBest(3) > maxDeflAllowed || ...
        changeBest(4) > maxDeflAllowed || ...
        changeBest(5) > maxDeflAllowed || ...
        changeBest(6) > maxDeflAllowed
    disp('howell_2 non è ok')
end


% Choice Blade Thickness in % of the chord (for both rotors)
thick1 = 0.1;      
thick2 = 0.1;

%Iteration on CAMBER ANGLE --> theta
i = 0;
maxiter = 100;
tol = 1e-3;
err = 0;

while err > tol && i < maxiter || i == 0

% epsilon = deltaBeta = theta + i - delta

% Link between CL and camber angle
% FOR NACA 65
% Cl = @(theta) theta / 25;
% FOR Circular arc profiles
Cl = @(theta) tand(theta/4) / 0.1103;

% Lift Coefficient evaluation
HUB.Cl1 = Cl(HUB.theta1);
MID.Cl1 = Cl(MID.theta1);
TIP.Cl1 = Cl(TIP.theta1);

HUB.Cl2 = Cl(HUB.theta2);
MID.Cl2 = Cl(MID.theta2);
TIP.Cl2 = Cl(TIP.theta2);

% Lieblein approach to evaluate the camber angle (through optimal incidence
% and deviation angle)
lieblein_rot1;
lieblein_rot2;

end

% Evaluate diffusion factor
checkLoading_rot1;
checkLoading_rot2;


% Choose Profiles for HUB, MID, TIP
profile = {'DCA', 'DCA', 'DCA',...
           'DCA', 'DCA', 'DCA'};

% Stagger angle evaluation
HUB.gamma1 = HUB.beta1 + HUB.i_opt1 + HUB.theta1 / 2;
MID.gamma1 = MID.beta1 + MID.i_opt1 + MID.theta1 / 2;
TIP.gamma1 = TIP.beta1 + TIP.i_opt1 + TIP.theta1 / 2;

HUB.gamma2 = HUB.beta3 + HUB.i_opt2 + HUB.theta2 / 2;
MID.gamma2 = MID.beta3 + MID.i_opt2 + MID.theta2 / 2;
TIP.gamma2 = TIP.beta3 + TIP.i_opt2 + TIP.theta2 / 2;


% Losses

% Profile Losses - Traupel Approach (with Mach and Re correction)
[HUB.Zprofile1, MID.Zprofile1, TIP.Zprofile1] = profileLosses_rot1('traupel',profile, HUB, MID, TIP);
[HUB.Zprofile2, MID.Zprofile2, TIP.Zprofile2] = profileLosses_rot2('traupel',profile, HUB, MID, TIP);

% Tip clearance / Leakage losses
deltaC = 2e-3; % tip leakage gap
[deltaEta_leakage1] = tipClearanceLoss(mdot,deltaC,b,TIP.c1,U1(2),MID.v1a,MID.leul1,TIP.beta1,TIP.beta2,eta1,S,rho1);
[deltaEta_leakage2] = tipClearanceLoss(mdot,deltaC,b,TIP.c2,U2(2),MID.v2a,MID.leul2,TIP.beta3,TIP.beta4,eta2,S,MID.rho2);

% Annulus Losses
% NO CORRELATION FOUND

% Secondary Losses / Endwall
t_TE1 = 0.01 * HUB.c1; %Trailing Edge Thickness
t_TE2 = 0.01 * HUB.c2;
[zetaEW1] = endwallLosses(HUB.Zprofile1,HUB.P2,HUB.Pt2,HUB.M2,HUB.deltaBeta1,rho1,HUB.rho2,HUB.v1,HUB.v2,HUB.alfa2,t_TE1,b,S);
[zetaEW2] = endwallLosses(HUB.Zprofile2,HUB.P4,HUB.Pt4,HUB.M4,HUB.deltaBeta2,HUB.rho2,HUB.rho4,HUB.v2,HUB.v4,HUB.alfa4,t_TE2,b,S);

% Disk Friction
delta = 2e-3; % gap between disks
[zetaD] = diskFrictionLosses(delta,Dhub,mdot,S,U1,U2,rho1,b,Cp,HUB,eta1,T1);

% Overall Losses Combination
% Yprofile --> ZetaP --> deltaEta_profile
[HUB.deltaEta_profile1] = deltaEta_calc('Z',[],HUB.Zprofile1,HUB.Mw2,P1,HUB.P2,T1,HUB.Beta1,Pt1,v1);
[MID.deltaEta_profile1] = deltaEta_calc('Z',[],MID.Zprofile1,MID.Mw2,P1,MID.P2,T1,MID.Beta1,Pt1,v1);
[TIP.deltaEta_profile1] = deltaEta_calc('Z',[],TIP.Zprofile1,TIP.Mw2,P1,TIP.P2,T1,TIP.Beta1,Pt1,v1);

[HUB.deltaEta_profile2] = deltaEta_calc('Z',[],HUB.Zprofile2,HUB.Mw4,HUB.P2,HUB.P4,HUB.T2,HUB.Beta2,HUB.Pt2,HUB.v2);
[MID.deltaEta_profile2] = deltaEta_calc('Z',[],MID.Zprofile2,MID.Mw4,MID.P2,MID.P4,MID.T2,MID.Beta2,MID.Pt2,MID.v2);
[TIP.deltaEta_profile2] = deltaEta_calc('Z',[],TIP.Zprofile2,TIP.Mw4,TIP.P2,TIP.P4,TIP.T2,TIP.Beta2,TIP.Pt2,TIP.v2);

% zetaEW --> deltaEta_endwall
[deltaEta_endwall1] = deltaEta_calc('Z',[],zetaEW1,HUB.M2,P1,HUB.P2,T1,HUB.Beta1,Pt1,HUB.v1);
[deltaEta_endwall2] = deltaEta_calc('Z',[],zetaEW2,HUB.M4,HUB.P2,HUB.P4,HUB.T2,HUB.Beta2,HUB.Pt2,HUB.v2);

% zetaD --> deltaEta_diskFriction
[deltaEta_diskFriction] = deltaEta_calc('Z',[],zetaD,[],P1,[],T1,BetaTot,Pt1,v1);

% Efficiencies estimated using losses correlations
eta1 = 1-((HUB.deltaEta_profile1 + MID.deltaEta_profile1 + TIP.deltaEta_profile1)/3 + deltaEta_endwall1 + deltaEta_leakage1);
eta2 = 1-((HUB.deltaEta_profile2 + MID.deltaEta_profile2 + TIP.deltaEta_profile2)/3 + deltaEta_endwall2 + deltaEta_leakage2);

% Total-to-Total Efficiency
deltaHtis = Cp * T1 * (BetaTot^GAMMA - 1) + (MID.v4^2 - MID.v1^2)/2;
etaTT = deltaHtis/leulTot;

% New Stage Efficiency
etaTOT= (( (MID.Beta1*MID.Beta2)^(GAMMA) - 1)*T1*eta1*eta2 ) / (   T1*eta2*(MID.Beta1^GAMMA-1) +  MID.T2*eta1*(MID.Beta2^GAMMA-1)   ) - deltaEta_diskFriction;

% Check on efficiency
err_efficiency = abs(etaTT - etaTOT)/etaTT;

i_eta = i_eta + 1;

end

% IGV

% Evaluate different IGV configurations for OFF Design requirements (reduced mass flow rate)
IGVoffDesign;

% Keep best configuration in terms of efficiency
IGVbest;

% Evaluate IGV characteristics during on Design operations
IGV;

% Plot Velocity Triangles
plotVelocities(HUBfp.v1a,HUBfp.v2a,HUBfp.v4a,HUBfp.v1,HUBfp.v2,HUBfp.v4,HUBfp.w1,HUBfp.w2,HUBfp.w3,HUBfp.w4,HUBfp.v1t,HUBfp.v2t,HUBfp.v4t,HUBfp.w1t,HUBfp.w2t,HUBfp.w3t,HUBfp.w4t,U1(1),U2(1),[1 0 0])
plotVelocities(MIDfp.v1a,MIDfp.v2a,MIDfp.v4a,MIDfp.v1,MIDfp.v2,MIDfp.v4,MIDfp.w1,MIDfp.w2,MIDfp.w3,MIDfp.w4,MIDfp.v1t,MIDfp.v2t,MIDfp.v4t,MIDfp.w1t,MIDfp.w2t,MIDfp.w3t,MIDfp.w4t,U1(2),U2(2),[0 1 0])
plotVelocities(TIPfp.v1a,TIPfp.v2a,TIPfp.v4a,TIPfp.v1,TIPfp.v2,TIPfp.v4,TIPfp.w1,TIPfp.w2,TIPfp.w3,TIPfp.w4,TIPfp.v1t,TIPfp.v2t,TIPfp.v4t,TIPfp.w1t,TIPfp.w2t,TIPfp.w3t,TIPfp.w4t,U1(3),U2(3),[0 0 1])

% Uncomment to plot 3D blades
%stampapalette

% Plot blades and writes .txt files containing blade points
bladePOINTS;

% Output data
% N. Blades 1 rotor
nBlades1
% N. Blades 2 rotor
nBlades2
% N. Blades IGV
nBladesIGV

% Mean Compression Ratio 1st rotor
Beta1=(MID.Beta1+HUB.Beta1+TIP.Beta1)/3
% Mean Compression Ratio 2nd rotor
Beta2=(MID.Beta2+HUB.Beta2+TIP.Beta2)/3
% Global Compression Ratio
Beta=(MID.Beta1*MID.Beta2+HUB.Beta1*HUB.Beta2+TIP.Beta1*TIP.Beta2)/3

% Efficiency 1st rotor
eta1
% Efficiency 2nd rotor
eta2
% Global Efficiency
etaTOT

% Global Compression Ratio - OFF DESIGN CONDITIONS
BetaoffDesign=(MIDfp.Beta1*MIDfp.Beta2+HUBfp.Beta1*HUBfp.Beta2+TIPfp.Beta1*TIPfp.Beta2)/3

% Efficiency - OFF DESIGN CONDITIONS
etaTOToffDesign = best_etaTOT

save results