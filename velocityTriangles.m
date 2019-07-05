function [station leulTot leul1 leul2]=velocityTriangles(v1, v1ax, v2ax, v4ax,S,P1, T1, U1, U2, etaTT,location,Dhis,work1,leul1,leul2,varargin)

    v1t=0;
    Cp=1004.69;
    R=287;
    gamma=1.4;
    mdot=100;

maxiter = 1000;    
    
v4=150; v2=240; 
v2old=1000; v4old=1000;

leulTot=Dhis/etaTT;

iterWork=0;     iterVelocity=[];     errorLeul=1000;

while errorLeul>1e-4
    
  iterWork=iterWork+1;
  leulTotOld=leulTot;
  
     if nargin<15
     leul1=leulTot*work1;
     leul2=leulTot-leul1;
     end
        
iterVelocity=[iterVelocity 0];   v2axold=1000;   v4old=1000;   v2old=1000;

while abs(v4-v4old)>1e-4 && abs(v2-v2old)>1e-4 && abs(v2ax-v2axold)>1e-4 && iterVelocity(end) < maxiter
    v2old=v2;
    v2axold=v2ax;
    v4old=v4;
    iterVelocity(end)=iterVelocity(end)+1;

    v2t=leul1/U1;
    v4t=leul2/U2+v2t;
    
    v2=sqrt(v2ax^2+v2t^2);
    v4=sqrt(v4ax^2+v4t^2);
    
    Dhis1=leul1*etaTT-(v2^2-v1^2)/2;
    Dhis2=leul2*etaTT-(v4^2-v2^2)/2;
    
    T2=T1+Dhis1/Cp/etaTT;
    T4=T2+Dhis2/Cp/etaTT;
    
    Beta1=(1+Dhis1/(Cp*T1))^(gamma/(gamma-1));
    Beta2=(1+Dhis2/(Cp*T2))^(gamma/(gamma-1));
    
    P2=Beta1*P1;
    P4=Beta2*P2;

    rho2=P2/(R*T2);
    rho4=P4/(R*T4);
    
    v2ax=mdot/(S*rho2);
    v4ax=mdot/(S*rho4);
end
    leulTot=(Dhis+(v4^2-v1^2)/2)/etaTT;
    
    errorLeul=abs(leulTot-leulTotOld);
    
    if nargin==15
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

w1=sqrt(w1t^2+v1ax^2);
w2=sqrt(w2t^2+v2ax^2);
w3=sqrt(w3t^2+v2ax^2);
w4=sqrt(w4t^2+v4ax^2);

beta1=atand(w1t/v1ax);
beta2=atand(w2t/v2ax);
beta3=atand(w3t/v2ax);
beta4=atand(w4t/v4ax);
Dbeta1=beta2-beta1;
Dbeta2=beta4-beta3;

alfa1=asind(v1t/v1);
alfa2=asind(v2t/v2);
alfa4=asind(v4t/v4);

plotVelocities(v1ax,v2ax,v4ax,v1,v2,v4,w1,w2,w3,w4,v1t,v2t,v4t,w1t,w2t,w3t,w4t,U1,U2,location)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
station.v1=v1;
station.v2=v2;
station.v4=v4;
station.v2ax=v2ax;
station.v4ax=v4ax;
station.P2=P2;
station.T2=T2;
station.rho2=rho2;
station.P4=P4;
station.T4=T4;
station.rho4=rho4;
station.v1t=v1t;
station.v2t=v2t;
station.v4t=v4t;
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
station.Dbeta1=Dbeta1;
station.Dbeta2=Dbeta2;
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

