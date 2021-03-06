function [best] = contraROTTO(param)
% DATI
mdot=100; %[kg/s]
Pt1=100000; %[Pa]
Tt1=300; %[K]
gamma=1.4; %Specific heat ratio
GAMMA=(gamma-1)/gamma;
R=287;
Cp=1004.69;
BetaTot=1.45; %compression ratio

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARAMETRI
Dtip=1;
eta1=0.76;
eta2=0.87;

work1=param(1); %lavoro percentuale sul primo rotore
rotRatio=param(2);   %se rotRatio cresce leul2 cresce ma M2rel cresce
n=param(3); %se n cresce leul1 cresce, Beta1 deve essere aumentato per non essere nullo all'hub
Mw1_tip=param(4);
Re1 = param(5);
Re2 = param(6);
alfa_mid = param(7); 
alfa_hub = param(8);
alfa_tip = param(9);

% alfa_hub = alfa_mid;
% alfa_tip = alfa_mid;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%
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


%% Outer Iteration on Efficiency
i_eta = 0;
maxiter_eta = 100;
tol = 1e-2;
err_efficiency = 0;

while err_efficiency > tol && i_eta < maxiter_eta || i_eta == 0
nBlades1 = 1;
nBlades2 = 1;
    
[MID leulTot leul1 leul2]=velocityTriangles(mdot,alfa_mid,v1,v1a_mid,v2a,v4a,S,P1,T1,U1(2),U2(2),eta1,eta2,[0 1 0],Dhis,work1);
if ~isreal(MID.v2t) 
    best = 100;
    return
end
[HUB]=velocityTriangles(mdot,alfa_hub,v1,v1a_hub,MID.v2a,MID.v4a,S,P1,T1,U1(1),U2(1),eta1,eta2,[1 0 0],leulTot,work1,leul1,leul2);
if ~isreal(HUB.v2t)
    best = 100;
    return
end
[TIP]=velocityTriangles(mdot,alfa_tip,v1,v1a_tip,MID.v2a,MID.v4a,S,P1,T1,U1(3),U2(3),eta1,eta2,[0 0 1],leulTot,work1,leul1,leul2);
if ~isreal(TIP.v2t)
    best = 100;
    return
end


if HUB.Beta1 < 1 || MID.Beta1 < 1 || TIP.Beta1 < 1 || HUB.Beta2 < 1 || MID.Beta2 < 1 || TIP.Beta2 < 1
    break
end
% check sui numeri di Mach
Mrel_tip1=TIP.w1/sqrt(gamma*R*T1);
Mrel_tip2=TIP.w3/sqrt(gamma*R*TIP.T2);


%Howell correlation per calcolare il
iterHowell = 0;
maxiter = 100;
tol = 1e-3;

sigmaLimit = 2.5;
sigma=0.625:0.01:sigmaLimit;


%while changeSolidity < 0.01 && changeSolidity > -0.01 || sigma > sigmaLimit && iterHowell < maxiter || iterHowell == 0

for i = 1 : length(sigma)
[deltaBetaOpt1(1)] = howellCorrelation(HUB.beta2,5e5,sigma(i));
[deltaBetaOpt1(2)] = howellCorrelation(MID.beta2,5e5,sigma(i));
[deltaBetaOpt1(3)] = howellCorrelation(TIP.beta2,5e5,sigma(i));
[deltaBetaOpt2(1)] = howellCorrelation(HUB.beta4,5e5,sigma(i));
[deltaBetaOpt2(2)] = howellCorrelation(MID.beta4,5e5,sigma(i));
[deltaBetaOpt2(3)] = howellCorrelation(TIP.beta4,5e5,sigma(i));

changeSolidity1(1)=deltaBetaOpt1(1)-HUB.deltaBeta1;  %se <1 aumentare solidity
changeSolidity1(2)=deltaBetaOpt1(2)-MID.deltaBeta1;  %se >1 diminuire solidity
changeSolidity1(3)=deltaBetaOpt1(3)-TIP.deltaBeta1;
changeSolidity2(1)=deltaBetaOpt2(1)+HUB.deltaBeta2;
changeSolidity2(2)=deltaBetaOpt2(2)+MID.deltaBeta2;
changeSolidity2(3)=deltaBetaOpt2(3)+TIP.deltaBeta2;
changeSolidity(:,i) = [changeSolidity1'; changeSolidity2'];

end

[changeBest,index] = min(abs(changeSolidity)');
maxDeflAllowed = 3;
if changeBest(1) > maxDeflAllowed || ...
        changeBest(2) > maxDeflAllowed || ...
        changeBest(3) > maxDeflAllowed || ...
        changeBest(4) > maxDeflAllowed || ...
        changeBest(5) > maxDeflAllowed || ...
        changeBest(6) > maxDeflAllowed
    disp('intervieni, howell non è ok')

    Mrel_tip2
    solidityBest = sigma(index);
    %return
else
solidityBest = sigma(index);
Mrel_tip2
end

%% New Part



mu = 1.81e-5;

%number of blades
[HUB, MID, TIP, nBlades1] = nBladesFUNC(HUB,MID,TIP,1,rho1,mu,Re1,solidityBest,Dhub,Dmid,Dtip);
[HUB, MID, TIP, nBlades2] = nBladesFUNC(HUB,MID,TIP,2,[],mu,Re2,solidityBest,Dhub,Dmid,Dtip);

nBlades1
nBlades2

if nBlades1 > 70 || nBlades2 > 70 || nBlades1 < 25 || nBlades2 < 25
    break
end

%NEW Howell

[deltaBetaOpt1(1)] = howellCorrelation(HUB.beta2,5e5,HUB.sigma1);
[deltaBetaOpt1(2)] = howellCorrelation(MID.beta2,5e5,MID.sigma1);
[deltaBetaOpt1(3)] = howellCorrelation(TIP.beta2,5e5,TIP.sigma1);
[deltaBetaOpt2(1)] = howellCorrelation(HUB.beta4,5e5,HUB.sigma2);
[deltaBetaOpt2(2)] = howellCorrelation(MID.beta4,5e5,MID.sigma2);
[deltaBetaOpt2(3)] = howellCorrelation(TIP.beta4,5e5,TIP.sigma2);

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
maxDeflAllowed = 3;
if changeBest(1) > maxDeflAllowed || ...
        changeBest(2) > maxDeflAllowed || ...
        changeBest(3) > maxDeflAllowed || ...
        changeBest(4) > maxDeflAllowed || ...
        changeBest(5) > maxDeflAllowed || ...
        changeBest(6) > maxDeflAllowed
    disp('intervieni, howell non è ok')
end


% CHOICE!!!!!!!!!!!!!!!!
thick1 = 0.1;      %thickness 10% of the chord
thick2 = 0.1;


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
Cl = @(theta) tand(theta/4) / 0.1103;

%theta is the camber angle
HUB.Cl1 = Cl(HUB.theta1);
MID.Cl1 = Cl(MID.theta1);
TIP.Cl1 = Cl(TIP.theta1);

HUB.Cl2 = Cl(HUB.theta2);
MID.Cl2 = Cl(MID.theta2);
TIP.Cl2 = Cl(TIP.theta2);

%Re = @(rho,u,d) rho * u * d / mu;
%
%Re_hub = Re((rho1+rho2_hub)/2,(w1_hub + w2_hub)/2,HUB.c1);
%Re_mid = Re((rho1+rho2_mid)/2,(w1_mid + w2_mid)/2,c);
%Re_tip = Re((rho1+rho2_tip)/2,(w1_tip + w2_tip)/2,c);
%keyboard
lieblein_rot1;
lieblein_rot2;

end

%keyboard
checkLoading_rot1;
checkLoading_rot2;


%% Choose Profiles for HUB, MID, TIP
profile = {'DCA', 'DCA', 'DCA',...
           'DCA', 'DCA', 'DCA'};

%stagger angle
HUB.gamma1 = HUB.beta1 + HUB.i_opt1 + HUB.theta1 / 2;
MID.gamma1 = MID.beta1 + MID.i_opt1 + MID.theta1 / 2;
TIP.gamma1 = TIP.beta1 + TIP.i_opt1 + TIP.theta1 / 2;

HUB.gamma2 = HUB.beta3 + HUB.i_opt2 + HUB.theta2 / 2;
MID.gamma2 = MID.beta3 + MID.i_opt2 + MID.theta2 / 2;
TIP.gamma2 = TIP.beta3 + TIP.i_opt2 + TIP.theta2 / 2;



%% Losses

etaTT_new = 1;
losses = 0;

%% Profile Losses
%profileLosses;

[HUB.Zprofile1, MID.Zprofile1, TIP.Zprofile1] = profileLosses_rot1('traupel',profile, HUB, MID, TIP);
[HUB.Zprofile2, MID.Zprofile2, TIP.Zprofile2] = profileLosses_rot2('traupel',profile, HUB, MID, TIP);

%% Tip clearance / Leakage losses
deltaC = 2e-3; % tip leakage gap
[deltaEta_leakage1] = tipClearanceLoss(mdot,deltaC,b,TIP.c1,U1(2),MID.v1a,MID.leul1,TIP.beta1,TIP.beta2,eta1,S,rho1);
[deltaEta_leakage2] = tipClearanceLoss(mdot,deltaC,b,TIP.c2,U2(2),MID.v2a,MID.leul2,TIP.beta3,TIP.beta4,eta2,S,MID.rho2);


%% Annulus Losses
% NO ANNULUS

%% Secondary Losses / Endwall
t_TE1 = 0.01 * HUB.c1;
t_TE2 = 0.01 * HUB.c2;
[zetaEW1] = endwallLosses(HUB.Zprofile1,HUB.P2,HUB.Pt2,HUB.M2,HUB.deltaBeta1,rho1,HUB.rho2,HUB.v1,HUB.v2,HUB.alfa2,t_TE1,b,S);
[zetaEW2] = endwallLosses(HUB.Zprofile2,HUB.P4,HUB.Pt4,HUB.M4,HUB.deltaBeta2,HUB.rho2,HUB.rho4,HUB.v2,HUB.v4,HUB.alfa4,t_TE2,b,S);

%% Disk Friction
delta = 2e-3; % gap between disks

[zetaD] = diskFrictionLosses(delta,Dhub,mdot,S,U1,U2,rho1,b,Cp,HUB,eta1,T1);

%% Overall Losses Combination (this fromula from Y to Zeta then to deltaEta)
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

[deltaEta_diskFriction] = deltaEta_calc('Z',[],zetaD,[],P1,[],T1,BetaTot,Pt1,v1);

eta1 = 1-((HUB.deltaEta_profile1 + MID.deltaEta_profile1 + TIP.deltaEta_profile1)/3 + deltaEta_endwall1 + deltaEta_leakage1);
eta2 = 1-((HUB.deltaEta_profile2 + MID.deltaEta_profile2 + TIP.deltaEta_profile2)/3 + deltaEta_endwall2 + deltaEta_leakage2);


%% Total-to-Total Efficiency
deltaHtis = Cp * T1 * (BetaTot^GAMMA - 1) + (MID.v4^2 - MID.v1^2)/2;
etaTT = deltaHtis/leulTot;

%% Check on efficiency
etaTOT= (( (MID.Beta1*MID.Beta2)^(GAMMA) - 1)*T1*eta1*eta2 ) / (   T1*eta2*(MID.Beta1^GAMMA-1) +  MID.T2*eta1*(MID.Beta2^GAMMA-1)   ) - deltaEta_diskFriction;
err_efficiency = abs(etaTT - etaTOT)/etaTT;

i_eta = i_eta + 1;

if eta1 < 0 || eta2 < 0
  break
end

end


if nBlades1 > 70 || nBlades2 > 70 || nBlades1 < 25 || nBlades2 < 25
    best = 10;
    return
    
elseif i_eta==maxiter_eta
    best = 9;
    return

elseif TIP.Mw3 > 0.88 || TIP.Mw1 > 0.88 || max(changeBest) > 5 ...
       || HUB.Mw3 > 0.88 || MID.Mw3 > 0.88
    best = 8;
    return

elseif eta1 < 0 || eta2 < 0 || eta1 > 1 || eta2 > 1 ...
        || HUB.Beta1 < 1 || MID.Beta1 < 1 || TIP.Beta1 < 1 ...
        || HUB.Beta2 < 1 || MID.Beta2 < 1 || TIP.Beta2 < 1
  best = 7;
  return

else    
    best = 1-etaTOT + 0.5*max(changeBest)+0.9*(HUB.Mw1+MID.Mw1+TIP.Mw1)+0.5*(HUB.Mw2+MID.Mw2);
    fprintf('eta1 = %f',eta1)
    fprintf('\neta2 = %f',eta2)
    fprintf('\netaTOT = %f',etaTOT)
end

end
