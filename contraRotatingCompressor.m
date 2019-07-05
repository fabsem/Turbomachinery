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
Mw1_tip=0.7;
etaTT=0.85;

work1=0.55; %lavoro percentuale sul primo rotore
rotRatio=1;   %se rotRatio cresce leul2 cresce ma M2rel cresce
n=3000; %se n cresce leul1 cresce, Beta1 deve essere aumentato per non essere nullo all'hub

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
Dhub=fsolve(Dhub,Dtip);
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

%Howell correlation per calcolare il
sigma1=[3 0.8 0.6];
sigma2=[3 2 1.1];
[DBetaOpt1(1)] = howellCorrelation(HUB.beta2,5e5,sigma1(1));
[DBetaOpt1(2)] = howellCorrelation(MID.beta2,5e5,sigma1(2));
[DBetaOpt1(3)] = howellCorrelation(TIP.beta2,5e5,sigma1(3));
[DBetaOpt2(1)] = howellCorrelation(HUB.beta4,5e5,sigma2(1));
[DBetaOpt2(2)] = howellCorrelation(MID.beta4,5e5,sigma2(2));
[DBetaOpt2(3)] = howellCorrelation(TIP.beta4,5e5,sigma2(3));

changeSolidity1(1)=DBetaOpt1(1)-HUB.Dbeta1;  %se <1 aumentare solidity
changeSolidity1(2)=DBetaOpt1(2)-MID.Dbeta1;  %se >1 diminuire solidity
changeSolidity1(3)=DBetaOpt1(3)-TIP.Dbeta1;
changeSolidity2(1)=DBetaOpt2(1)+HUB.Dbeta2;
changeSolidity2(2)=DBetaOpt2(2)+MID.Dbeta2;
changeSolidity2(3)=DBetaOpt2(3)+TIP.Dbeta2;

%     Dh1is_tip=mid.DhTT1is-(tip.v2^2-tip.v1^2)/2;
%     Dh1is_hub=mid.DhTT1is-(hub.v2^2-hub.v1^2)/2;
%     Dh2is_tip=mid.DhTT2is-(tip.v4^2-tip.v2^2)/2;
%     Dh2is_hub=mid.DhTT2is-(hub.v4^2-hub.v2^2)/2;
%         
% Dh1is=[Dh1is_hub   Cp*T1*(Beta1(2)^((gamma-1)/gamma)-1)   Dh1is_tip];  %hub-mid-tip
% Dh2is=[Dh2is_hub   Cp*T2(2)*(Beta2(2)^((gamma-1)/gamma)-1)   Dh2is_tip];
% 
%     clearvars Dh1is_tip Dh1is_hub Dh2is_tip Dh2is_hub
% 
% T2=[T1+Dh1is(1)/Cp/etaTT     T1+Dh1is(2)/Cp/etaTT     T1+Dh1is(3)/Cp/etaTT] ; %hub-mid-tip
% T4=[T2(1)+Dh2is(1)/Cp/etaTT   T2(2)+Dh2is(2)/Cp/etaTT   T2(3)+Dh2is(3)/Cp/etaTT];
% 
% 
% Beta1=[(1+Dh1is(1)/(Cp*T1))^(gamma/(gamma-1))     Beta1(2)     (1+Dh1is(3)/(Cp*T1))^(gamma/(gamma-1))];
% Beta2=[(1+Dh2is(1)/(Cp*T2(1)))^(gamma/(gamma-1))     Beta2(2)     (1+Dh2is(3)/(Cp*T2(3)))^(gamma/(gamma-1))];
% P2=[Beta1(1)*P1    Beta1(2)*P1     Beta1(3)*P1];
% P4=[Beta2(1)*P2(1)    Beta2(2)*P2(2)     Beta2(3)*P2(3)];
% 
% rho2=[P2(1)/(R*T2(1))    P2(2)/(R*T2(2))    P2(3)/(R*T2(3))];
% rho4=[P4(1)/(R*T4(1))    P4(2)/(R*T4(2))    P4(3)/(R*T4(3))];
% 
% 
% v2ax=[mdot/(S*rho2(1))    mdot/(S*rho2(2))    mdot/(S*rho2(3))];
% v4ax=[mdot/(S*rho4(1))    mdot/(S*rho4(2))    mdot/(S*rho4(3))];
% 


%plotVelocities(v1ax,v1,v2,v4,w1,w2,w3,w4,v1t,v2t,v4t,w1t,w2t,w3t,w4t,U1_mid,U2_mid)

%leul_check=0.5*(v2^2-v1^2)-0.5*(w2^2-w1^2)+0.5*(v4^2-v2^2)-0.5*(w4^2-w3^2);




