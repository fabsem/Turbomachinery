% IGV

mdot=90;

MIDfp.v1a=mdot/(S*rho1);
MIDfp.w1a=MIDfp.v1a;
MIDfp.v1 = MIDfp.v1a;
alfa1 = 0;

%Guess Beta_FP and try to match the reduced mass flow rate
Beta_new=1.4;
i_mdotfp = 0;
err_deltaBeta = 0;
maxiter = 100000;
tol = 1e-2;
MIDfp.v2a=MIDfp.v1a;
MIDfp.v4a=MIDfp.v1a;
mdot_iter=mdot;

while err_deltaBeta > tol && i_mdotfp < maxiter || i_mdotfp == 0
BetaTot = Beta_new;

Dhis=Cp*T1*(BetaTot^((gamma-1)/gamma)-1);


[MIDfp leulTotfp leul1fp leul2fp]=velocityTriangles(mdot_iter,alfa1,MIDfp.v1, MIDfp.v1a, MIDfp.v2a, MIDfp.v4a,S,P1,T1,U1(2),U2(2), etaTT,[0 1 0],Dhis,work1);


err_deltaBeta = abs(MID.deltaBeta1-MIDfp.deltaBeta1);

i_mdotfp = i_mdotfp + 1;
Beta_new = BetaTot + 0.0001;
if Beta_new > 2
  disp('Non c''è soluzione')
  break
 end
end

BetaTot