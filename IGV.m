% IGV

mdot=90;

v1a=mdot/(S*rho1);
w1a=v1a;

HUBfp.w1t=tand(HUB.beta1)*w1a;
MIDfp.w1t=tand(MID.beta1)*w1a;
TIPfp.w1t=tand(TIP.beta1)*w1a;

HUBfp.v1t=HUBfp.w1t+U1;
MIDfp.v1t=MIDfp.w1t+U1;
TIPfp.v1t=TIPfp.w1t+U1;

HUBfp.alfa0=atan(HUBfp.v1t/v1a);
MIDfp.alfa0=atan(MIDfp.v1t/v1a);
TIPfp.alfa0=atan(TIPfp.v1t/v1a);
