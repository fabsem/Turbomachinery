%PATH
addpath('liebleinCorrelations/')
addpath('howellCorrelations/')
addpath('losses/')
addpath('IGV/')
% DATI
clear all; close all; clc
mdot=100; %[kg/s]
Pt1=100000; %[Pa]
Tt1=300; %[K]
gamma=1.4; %Specific heat ratio
R=287;
Cp=1004.69;
BetaTot=1.45; %compression ratio

%Inlet guide vane is REQUIRED
%Free to choose rpm
%DESIGN WITH HIGHEST ACHIEVABLE EFFICIENCY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARAMETRI
Dtip=1;
etaTT=0.85;


Mw1_tip=0.773;%0.7823;
work1=0.4766;%0.5;%0.3;%0.3330;%0.3317;%0.3387;%0.339;%0.3177;%0.3; %lavoro percentuale sul primo rotore
rotRatio=0.78;%1.2455;%1.2574;%1.5;%1.3;   %se rotRatio cresce leul2 cresce ma M2rel cresce
n=3000;%2870;%2800;%2848;%3000; %se n cresce leul1 cresce, Beta1 deve essere aumentato per non essere nullo all'hub

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
    w1t_tip=-U1_tip;
    v1a=v1;
    w1_tip=sqrt(v1a^2+w1t_tip^2);
    Mrel_tip_check=w1_tip/sqrt(gamma*R*T1);
    residual=abs(Mw1_tip-Mrel_tip_check);
end
%%
clearvars U1_tip residual w1_tip w1t_tip

P1=Pt1*(T1/Tt1)^(gamma/(gamma-1));
rho1=P1/(R*T1);
S=mdot/(v1a*rho1);
Dhub=@(Dhub) S-pi/4*(Dtip^2-Dhub^2);
Dhub=fzero(Dhub,Dtip);
b=(Dtip-Dhub)/2;
Dmid=(Dtip+Dhub)/2;

U1=[omega*Dhub/2 omega*Dmid/2 omega*Dtip/2];
U2=[-U1(1)*rotRatio -U1(2)*rotRatio -U1(3)*rotRatio]; %solo se b1=b2

Dhis=Cp*T1*(BetaTot^((gamma-1)/gamma)-1);
v2a=v1a; v4a=v1a; alpha=0;

[MID leulTot leul1 leul2]=velocityTriangles(mdot,alpha,v1,v1a,v2a,v4a,S,P1,T1,U1(2),U2(2),etaTT,[0 1 0],Dhis,work1);
[HUB]=velocityTriangles(mdot,alpha,v1,v1a,MID.v2a,MID.v4a,S,P1,T1,U1(1),U2(1),etaTT,[1 0 0],leulTot,work1,leul1,leul2);
[TIP]=velocityTriangles(mdot,alpha,v1,v1a,MID.v2a,MID.v4a,S,P1,T1,U1(3),U2(3),etaTT,[0 0 1],leulTot,work1,leul1,leul2);

% check sui numeri di Mach
Mrel_tip1=TIP.w1/sqrt(gamma*R*T1);
Mrel_tip2=TIP.w3/sqrt(gamma*R*TIP.T2);

close all
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
maxDeflAllowed = 2;
if changeBest(1) > maxDeflAllowed || ...
        changeBest(2) > maxDeflAllowed || ...
        changeBest(3) > maxDeflAllowed || ...
        changeBest(4) > maxDeflAllowed || ...
        changeBest(5) > maxDeflAllowed || ...
        changeBest(6) > maxDeflAllowed
    disp('intervieni, howell non è ok')

    Mrel_tip2
    solidityBest = sigma(index)
    return
else
solidityBest = sigma(index)
Mrel_tip2
end

%% New Part


Re1 = 6e5;
Re2 = 1.5e6;
mu = 1.81e-5;

%number of blades
[HUB, MID, TIP, nBlades1] = nBladesFUNC(HUB,MID,TIP,1,rho1,mu,Re1,solidityBest,Dhub,Dmid,Dtip);
[HUB, MID, TIP, nBlades2] = nBladesFUNC(HUB,MID,TIP,2,[],mu,Re2,solidityBest,Dhub,Dmid,Dtip);

nBlades1
nBlades2


%NEW Howell
for i = 1 : length(sigma)
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
changeSolidity(:,i) = [changeSolidity1'; changeSolidity2'];

end

[changeBest,~] = min(abs(changeSolidity)')
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
%chord1 = 0.13;       %chord [m] constant over the span
%chord2 = 0.13;

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
Cl = @(theta) theta / 25;

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

lieblein_rot1;
lieblein_rot2;

end

%keyboard
checkLoading_rot1;
checkLoading_rot2;


%% Losses
%losses;

%% IGV
IGV;
