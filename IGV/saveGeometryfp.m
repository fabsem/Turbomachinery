
HUBfp.c1=HUB.c1;
MIDfp.c1=MID.c1;
TIPfp.c1=TIP.c1;

HUBfp.c2=HUB.c2;
MIDfp.c2=MID.c2;
TIPfp.c2=TIP.c2;

HUBfp.sigma1=HUB.sigma1;
MIDfp.sigma1=MID.sigma1;
TIPfp.sigma1=TIP.sigma1;

HUBfp.sigma2=HUB.sigma2;
MIDfp.sigma2=MID.sigma2;
TIPfp.sigma2=TIP.sigma2;

HUBfp.theta_c1=HUB.theta_c1;
MIDfp.theta_c1=MID.theta_c1;
TIPfp.theta_c1=TIP.theta_c1;

HUBfp.theta_c2=HUB.theta_c2;
MIDfp.theta_c2=MID.theta_c2;
TIPfp.theta_c2=TIP.theta_c2;

HUBfp.M2 = HUBfp.v2/sqrt(gamma*R*HUBfp.T2);
MIDfp.M2 = MIDfp.v2/sqrt(gamma*R*MIDfp.T2);
TIPfp.M2 = TIPfp.v2/sqrt(gamma*R*TIPfp.T2);

HUBfp.M4 = HUBfp.v4/sqrt(gamma*R*HUBfp.T4);
MIDfp.M4 = MIDfp.v4/sqrt(gamma*R*MIDfp.T4);
TIPfp.M4 = TIPfp.v4/sqrt(gamma*R*TIPfp.T4);

HUBfp.Pt2 = HUBfp.P2*(1+(gamma-1)/2*HUBfp.M2^2)^(gamma/(gamma-1));
MIDfp.Pt2 = MIDfp.P2*(1+(gamma-1)/2*MIDfp.M2^2)^(gamma/(gamma-1));
TIPfp.Pt2 = TIPfp.P2*(1+(gamma-1)/2*TIPfp.M2^2)^(gamma/(gamma-1));

HUBfp.Pt4 = HUBfp.P4*(1+(gamma-1)/2*HUBfp.M4^2)^(gamma/(gamma-1));
MIDfp.Pt4 = MIDfp.P4*(1+(gamma-1)/2*MIDfp.M4^2)^(gamma/(gamma-1));
TIPfp.Pt4 = TIPfp.P4*(1+(gamma-1)/2*TIPfp.M4^2)^(gamma/(gamma-1));

HUBfp.w1 = sqrt(HUBfp.w1t^2+HUBfp.v1a^2);
MIDfp.w1 = sqrt(MIDfp.w1t^2+MIDfp.v1a^2);
TIPfp.w1 = sqrt(TIPfp.w1t^2+TIPfp.v1a^2);

HUBfp.w2 = sqrt(HUBfp.w2t^2+HUBfp.v2a^2);
MIDfp.w2 = sqrt(MIDfp.w2t^2+MIDfp.v2a^2);
TIPfp.w2 = sqrt(TIPfp.w2t^2+TIPfp.v2a^2);

HUBfp.w3 = sqrt(HUBfp.w3t^2+HUBfp.v2a^2);
MIDfp.w3 = sqrt(MIDfp.w3t^2+MIDfp.v2a^2);
TIPfp.w3 = sqrt(TIPfp.w3t^2+TIPfp.v2a^2);

HUBfp.w4 = sqrt(HUBfp.w4t^2+HUBfp.v4a^2);
MIDfp.w4 = sqrt(MIDfp.w4t^2+MIDfp.v4a^2);
TIPfp.w4 = sqrt(TIPfp.w4t^2+TIPfp.v4a^2);


HUBfp.Mw1 = HUBfp.w1/sqrt(gamma*R*T1);
MIDfp.Mw1 = MIDfp.w1/sqrt(gamma*R*T1);
TIPfp.Mw1 = TIPfp.w1/sqrt(gamma*R*T1);

HUBfp.Mw2 = HUBfp.w2/sqrt(gamma*R*HUBfp.T2);
MIDfp.Mw2 = MIDfp.w2/sqrt(gamma*R*MIDfp.T2);
TIPfp.Mw2 = TIPfp.w2/sqrt(gamma*R*TIPfp.T2);

HUBfp.Mw3 = HUBfp.w3/sqrt(gamma*R*HUBfp.T2);
MIDfp.Mw3 = MIDfp.w3/sqrt(gamma*R*MIDfp.T2);
TIPfp.Mw3 = TIPfp.w3/sqrt(gamma*R*TIPfp.T2);

HUBfp.Mw4 = HUBfp.w4/sqrt(gamma*R*HUBfp.T4);
MIDfp.Mw4 = MIDfp.w4/sqrt(gamma*R*MIDfp.T4);
TIPfp.Mw4 = TIPfp.w4/sqrt(gamma*R*TIPfp.T4);

HUBfp.Re1 = rho1 * HUBfp.w1 * HUBfp.c1 / mu;
MIDfp.Re1 = rho1 * MIDfp.w1 * MIDfp.c1 / mu;
TIPfp.Re1 = rho1 * TIPfp.w1 * TIPfp.c1 / mu;

HUBfp.Re2 = HUBfp.rho2 * HUBfp.w3 * HUBfp.c2 / mu;
MIDfp.Re2 = MIDfp.rho2 * MIDfp.w3 * MIDfp.c2 / mu;
TIPfp.Re2 = TIPfp.rho2 * TIPfp.w3 * TIPfp.c2 / mu;
