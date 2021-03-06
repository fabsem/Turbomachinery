function [station leulTot leul1 leul2]=velocityTriangles(mdot,alfa1,v1, v1a, v2a, v4a,S,P1, T1, U1, U2, eta1,eta2,location,Dhis,work1,leul1,leul2,varargin)

    v1t=v1 * sind(alfa1);
    Cp=1004.69;
    R=287;
    gamma=1.4;
    GAMMA=(gamma-1)/gamma;
maxiter = 1000;

v4=150; v2=240;
v2old=1000; v4old=1000;

leulTot=Dhis/eta1;

iterWork=0;     iterVelocity=[];     errorLeul=1000;

while errorLeul>1e-4 && iterWork < maxiter

  iterWork=iterWork+1;
  leulTotOld=leulTot;

     if nargin<18
     leul1=leulTot*work1;
     leul2=leulTot-leul1;
     end

iterVelocity=[iterVelocity 0];   v2axold=1000;   v4old=1000;   v2old=1000;

while abs(v4-v4old)>1e-4 && abs(v2-v2old)>1e-4 && abs(v2a-v2axold)>1e-4 && iterVelocity(end) < maxiter
    v2old=v2;
    v2axold=v2a;
    v4old=v4;
    iterVelocity(end)=iterVelocity(end)+1;

    v2t=leul1/U1;
    v4t=leul2/U2+v2t;

    v2=sqrt(v2a^2+v2t^2);
    v4=sqrt(v4a^2+v4t^2);

    Dhis1=leul1*eta1-(v2^2-v1^2)/2;
    Dhis2=leul2*eta2-(v4^2-v2^2)/2;

    T2=T1+Dhis1/Cp/eta1;
    T4=T2+Dhis2/Cp/eta2;

    Beta1=(1+Dhis1/(Cp*T1))^(gamma/(gamma-1));
    Beta2=(1+Dhis2/(Cp*T2))^(gamma/(gamma-1));

    P2=Beta1*P1;
    P4=Beta2*P2;

    rho2=P2/(R*T2);
    rho4=P4/(R*T4);

    v2a=mdot/(S*rho2);
    v4a=mdot/(S*rho4);
end
    etaTOT= (( (Beta1*Beta2)^(GAMMA) - 1)*T1*eta1*eta2 ) / (   T1*eta2*(Beta1^GAMMA-1) +  T2*eta1*(Beta2^GAMMA-1)   );
    leulTot=(Dhis+(v4^2-v1^2)/2)/etaTOT;

    errorLeul=abs(leulTot-leulTotOld);

    if nargin==18
       errorLeul=0;
    end
end
if ~isreal(v2t)
    station.v2t=v2t;
    return
else
w1t=v1t-U1;
w2t=v2t-U1;
w3t=v2t-U2;
w4t=v4t-U2;

w1=sqrt(w1t^2+v1a^2);
w2=sqrt(w2t^2+v2a^2);
w3=sqrt(w3t^2+v2a^2);
w4=sqrt(w4t^2+v4a^2);

beta1=atand(w1t/v1a);
beta2=atand(w2t/v2a);
beta3=atand(w3t/v2a);
beta4=atand(w4t/v4a);
deltaBeta1=beta2-beta1;
deltaBeta2=beta4-beta3;

alfa1=asind(v1t/v1);
alfa2=asind(v2t/v2);
alfa4=asind(v4t/v4);

%plotVelocities(v1ax,v2ax,v4ax,v1,v2,v4,w1,w2,w3,w4,v1t,v2t,v4t,w1t,w2t,w3t,w4t,U1,U2,location)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
station.v1=v1;
station.v2=v2;
station.v4=v4;
station.v1a=v1a;
station.v2a=v2a;
station.v4a=v4a;
station.P2=P2;
station.T2=T2;
station.rho2=rho2;
station.P4=P4;
station.T4=T4;
station.rho4=rho4;
station.v1t=v1t;
station.v2t=v2t;
station.v4t=v4t;
station.Mw1=w1/sqrt(gamma*R*T1);
station.Mw2=w2/sqrt(gamma*R*T2);
station.Mw3=w3/sqrt(gamma*R*T2);
station.Mw4=w4/sqrt(gamma*R*T4);
station.M1=v1/sqrt(gamma*R*T1);
station.M2=v2/sqrt(gamma*R*T2);
station.M4=v4/sqrt(gamma*R*T4);
station.Pt2=P2*(1+(gamma-1)/2*station.M2^2)^(gamma/(gamma-1));
station.Pt4=P4*(1+(gamma-1)/2*station.M4^2)^(gamma/(gamma-1));
station.w1=w1;
station.w2=w2;
station.w3=w3;
station.w4=w4;
station.w1t=w1t;
station.w2t=w2t;
station.w3t=w3t;
station.w4t=w4t;
station.alfa1=alfa1;
station.alfa2=alfa2;
station.alfa4=alfa4;
station.beta1=beta1;
station.beta2=beta2;
station.beta3=beta3;
station.beta4=beta4;
station.deltaBeta1=deltaBeta1;
station.deltaBeta2=deltaBeta2;
station.leul1=leul1;
station.leul2=leul2;
station.Beta1=Beta1;
station.Beta2=Beta2;
station.Dhis1=Dhis1;
station.Dhis2=Dhis2;
station.iterWork=iterWork;
station.iterVelocity=iterVelocity;
end
end
