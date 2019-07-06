% IGV

mdot=90;

v1a=mdot/(S*rho1);
w1a=v1a;

HUBfp.w1t=tand(HUB.beta1)*w1a;
MIDfp.w1t=tand(MID.beta1)*w1a;
TIPfp.w1t=tand(TIP.beta1)*w1a;

HUBfp.v1t=HUBfp.w1t+U1(1);
MIDfp.v1t=MIDfp.w1t+U1(2);
TIPfp.v1t=TIPfp.w1t+U1(3);

HUBfp.v1=sqrt(HUBfp.v1t^2+v1a^2);
MIDfp.v1=sqrt(MIDfp.v1t^2+v1a^2);
TIPfp.v1=sqrt(TIPfp.v1t^2+v1a^2);

HUBfp.alfa1=atand(HUBfp.v1t/v1a);
MIDfp.alfa1=atand(MIDfp.v1t/v1a);
TIPfp.alfa1=atand(TIPfp.v1t/v1a);

%Guess Beta_FP and try to match the reduced mass flow rate
Beta_new=1.1;
i_mdotfp = 0;
err_deltaBeta = 0;
maxiter = 100000;
tol = 1e-4;
MIDfp.v2ax=v1a;
MIDfp.v4ax=v1a;
mdot_iter=mdot;

while err_deltaBeta > tol && i_mdotfp < maxiter || i_mdotfp == 0
BetaTot = Beta_new

Dhis=Cp*T1*(BetaTot^((gamma-1)/gamma)-1);


[MIDfp leulTotfp leul1fp leul2fp]=velocityTriangles(mdot_iter,MIDfp.alfa1,MIDfp.v1, v1a, MIDfp.v2ax, MIDfp.v4ax,S,P1, T1,U1(2),U2(2), etaTT,[0 1 0],Dhis,work1);


err_deltaBeta = abs(MID.deltaBeta1-MIDfp.deltaBeta1)/MID.deltaBeta1;

i_mdotfp = i_mdotfp + 1
Beta_new = BetaTot + 0.0001;
if Beta_new > 2
  disp('Non c''è soluzione')
  return
end
end
