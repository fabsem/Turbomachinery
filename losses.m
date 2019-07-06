%% Losses Rotor 1
%Profile losses

%Annulus losses
%HOWELL
deltaCDa = howell_AnnulusLosses(MID.s1,b);
%Vavra
deltaCDa = vavra_AnnulusLosses(MID.c1,b);

%Secondary Flow losses
%HOWELL
deltaCDs = howell_SecondaryLosses([],MID.Cl1);
%Vavra
deltaCDs = vavra_SecondaryLosses(MID.c1,b,MID.Cl1);

%Leakage/Recirculation losses
%Laksminaraiana
delta=0.025*b;
deltaCDdelta = laksminaraiana_LeakageLosses(delta,b,TIP.Cl1);
%Vavra
deltaCDdelta = vavra_LeakageLosses(delta,b,MID.c1,MID.s1,MID.alfa2,TIP.Cl1);

%% Losses Rotor 2
%Profile losses

%Annulus losses
%HOWELL
deltaCDa = howell_AnnulusLosses(MID.s2,b);
%Vavra
deltaCDa = vavra_AnnulusLosses(MID.c2,b);

%Secondary Flow losses
%HOWELL
deltaCDs = howell_SecondaryLosses([],MID.Cl2);
%Vavra
deltaCDs = vavra_SecondaryLosses(MID.c2,b,MID.Cl2);

%Leakage/Recirculation losses
%Laksminaraiana
delta=0.025*b;
deltaCDdelta = laksminaraiana_LeakageLosses(delta,b,TIP.Cl2);
%Vavra
deltaCDdelta = vavra_LeakageLosses(delta,b,MID.c2,MID.s2,MID.alfa4,TIP.Cl2);

%Disk Friction losses
