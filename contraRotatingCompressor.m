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
    v1ax=v1;
    w1_tip=sqrt(v1ax^2+w1t_tip^2);
    Mrel_tip_check=w1_tip/sqrt(gamma*R*T1);
    residual=abs(Mw1_tip-Mrel_tip_check);
end    
%%
clearvars U1_tip residual w1_tip w1t_tip

P1=Pt1*(T1/Tt1)^(gamma/(gamma-1));
rho1=P1/(R*T1);
S=mdot/(v1ax*rho1);
Dhub=@(Dhub) S-pi/4*(Dtip^2-Dhub^2);
Dhub=fzero(Dhub,Dtip);
b=(Dtip-Dhub)/2;
Dmid=(Dtip+Dhub)/2;

U1=[omega*Dhub/2 omega*Dmid/2 omega*Dtip/2];
U2=[-U1(1)*rotRatio -U1(2)*rotRatio -U1(3)*rotRatio]; %solo se b1=b2

Dhis=Cp*T1*(BetaTot^((gamma-1)/gamma)-1);
v2ax=v1ax; v4ax=v1ax;

[MID leulTot leul1 leul2]=velocityTriangles(v1,v1ax,v2ax,v4ax,S,P1,T1,U1(2),U2(2),etaTT,[0 1 0],Dhis,work1);
[HUB]=velocityTriangles(v1,v1ax,MID.v2ax,MID.v4ax,S,P1,T1,U1(1),U2(1),etaTT,[1 0 0],leulTot,work1,leul1,leul2);
[TIP]=velocityTriangles(v1,v1ax,MID.v2ax,MID.v4ax,S,P1,T1,U1(3),U2(3),etaTT,[0 0 1],leulTot,work1,leul1,leul2);

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
[DBetaOpt1(1)] = howellCorrelation(HUB.beta2,5e5,sigma(i));
[DBetaOpt1(2)] = howellCorrelation(MID.beta2,5e5,sigma(i));
[DBetaOpt1(3)] = howellCorrelation(TIP.beta2,5e5,sigma(i));
[DBetaOpt2(1)] = howellCorrelation(HUB.beta4,5e5,sigma(i));
[DBetaOpt2(2)] = howellCorrelation(MID.beta4,5e5,sigma(i));
[DBetaOpt2(3)] = howellCorrelation(TIP.beta4,5e5,sigma(i));

changeSolidity1(1)=DBetaOpt1(1)-HUB.Dbeta1;  %se <1 aumentare solidity
changeSolidity1(2)=DBetaOpt1(2)-MID.Dbeta1;  %se >1 diminuire solidity
changeSolidity1(3)=DBetaOpt1(3)-TIP.Dbeta1;
changeSolidity2(1)=DBetaOpt2(1)+HUB.Dbeta2;
changeSolidity2(2)=DBetaOpt2(2)+MID.Dbeta2;
changeSolidity2(3)=DBetaOpt2(3)+TIP.Dbeta2;
changeSolidity(:,i) = [changeSolidity1'; changeSolidity2'];

end

[changeBest,index] = min(abs(changeSolidity)')

if changeBest(1) > 1 || changeBest(2) > 1 || changeBest(3) > 1 || changeBest(4) > 1 || changeBest(5) > 1 || changeBest(6) > 1
    disp('intervieni, howell non è ok')
    
    Mrel_tip2
    solidityBest = sigma(index)
    return 
else
solidityBest = sigma(index)
Mrel_tip2
end
%iterHowell = iterHowell + 1;



%% New Part
Re = 3e5;
mu = 1.81e-5;

c_hub = Re * mu / (rho1 * w1_hub);
c_mid = Re * mu / (rho1 * w1_mid);
c_tip = Re * mu / (rho1 * w1_tip);

%if deltaBeta_optimal --> SOMETHING TO DO ON THE SOLIDITY?

%% Skip Howell, Pass to Leiblein
% Assumptions:
% etaStage = etaS = 0.85
% eta_rotor = etaStage
%
% Choose:
th = 0.08;      %thickness 8% of the chord
c = 0.13;       %chord [m] constant over the span

%number of blades
sigma_mid = 1;
s = c / sigma_mid;

n_blades = pi * D_mid / s;
n_blades = floor(n_blades);
if rem(n_blades,2)
    n_blades = n_blades + 1;
end

s_mid = pi * D_mid / n_blades;
s_tip = pi * D_tip / n_blades;
s_hub = pi * D_hub / n_blades;

sigma_mid = c / s_mid; sigma_mid = 1; %imposed!!!!!!
sigma_tip = c / s_tip; 
sigma_hub = c / s_hub;

%% Camber angle --> need iteration
% epsilon = deltaBeta = theta + i - delta
% DO GRAPHS FOR DIFFERENT ANGLES!!!
epsilon_hub = deltaBeta_hub;
epsilon_mid = deltaBeta_mid;
epsilon_tip = deltaBeta_tip;

%guess
theta_hub = epsilon_hub;%16;
theta_mid = epsilon_mid;%35;
theta_tip = epsilon_tip;%52;2;


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

Cl_hub = Cl(theta_hub);
Cl_mid = Cl(theta_mid);
Cl_tip = Cl(theta_tip);

Re = @(rho,u,d) rho * u * d / mu;

Re_hub = Re((rho1+rho2_hub)/2,(w1_hub + w2_hub)/2,c);
Re_mid = Re((rho1+rho2_mid)/2,(w1_mid + w2_mid)/2,c);
Re_tip = Re((rho1+rho2_tip)/2,(w1_tip + w2_tip)/2,c);


lieblein;

end


checkLoading;