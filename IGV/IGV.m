% IGV

mdot=90;

MIDfp.v1a=mdot/(S*rho1);
MIDfp.w1a=MIDfp.v1a;

HUBfp.w1t=tand(HUB.beta1)*MIDfp.w1a;
MIDfp.w1t=tand(MID.beta1)*MIDfp.w1a;
TIPfp.w1t=tand(TIP.beta1)*MIDfp.w1a;

HUBfp.v1t=HUBfp.w1t+U1(1);
MIDfp.v1t=MIDfp.w1t+U1(2);
TIPfp.v1t=TIPfp.w1t+U1(3);

HUBfp.v1=sqrt(HUBfp.v1t^2+MIDfp.v1a^2);
MIDfp.v1=sqrt(MIDfp.v1t^2+MIDfp.v1a^2);
TIPfp.v1=sqrt(TIPfp.v1t^2+MIDfp.v1a^2);

HUBfp.alfa1=atand(HUBfp.v1t/MIDfp.v1a);
alfa1=atand(MIDfp.v1t/MIDfp.v1a);
TIPfp.alfa1=atand(TIPfp.v1t/MIDfp.v1a);

%Guess Beta_FP and try to match the reduced mass flow rate
Beta_new=1.3;
i_mdotfp = 0;
err_deltaBeta = 0;
maxiter = 100000;
tol = 1e-1;
MIDfp.v2a=MIDfp.v1a;
MIDfp.v4a=MIDfp.v1a;
mdot_iter=mdot;

while err_deltaBeta > tol && i_mdotfp < maxiter || i_mdotfp == 0
BetaTot = Beta_new

Dhis=Cp*T1*(BetaTot^((gamma-1)/gamma)-1);


[MIDfp leulTotfp leul1fp leul2fp]=velocityTriangles(mdot_iter,alfa1,MIDfp.v1, MIDfp.v1a, MIDfp.v2a, MIDfp.v4a,S,P1,T1,U1(2),U2(2), etaTT,[0 1 0],Dhis,work1);
keyboard

err_deltaBeta = abs(MID.deltaBeta1-MIDfp.deltaBeta1);

i_mdotfp = i_mdotfp + 1
Beta_new = BetaTot + 0.0001;
if Beta_new > 2
  disp('Non c''è soluzione')
  return
 end
end
