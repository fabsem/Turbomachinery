% IGVonDesign

%% Geometrical Angle
% i = i0 + n*theta
% i0 = Kth * Kshape * i0.10

% DO GRAPHS --> [i010] = liblein_incidence010(beta1,sigma)
i010_tip = lieblein_incidence010(0,TIPfp.sigma_IGV);%4;
i010_mid = lieblein_incidence010(0,MIDfp.sigma_IGV);%4;
i010_hub = lieblein_incidence010(0,HUBfp.sigma_IGV);%4.2;

%correction coefficient for the shape
Ki_shape = lieblein_icorrectionShape('naca65'); %1; %for NACA 65

%correction for the thickness
% do graph --> [K_thick] = liblein_ithick(t,c)
%[Ki_thick] = lieblein_icorrectionThick(t_c)
Ki_thick_tip = lieblein_icorrectionThick(thick_IGV); %0.9;
Ki_thick_mid = lieblein_icorrectionThick(thick_IGV); %0.9;
Ki_thick_hub = lieblein_icorrectionThick(thick_IGV); %0.9;

% n coefficient
% do graph --> [n_lieblein] = lieblein_ncoeff(beta1,sigma)
n_lieblein_tip = lieblein_ncoeff(0,TIPfp.sigma_IGV); %-0.26;
n_lieblein_mid = lieblein_ncoeff(0,MIDfp.sigma_IGV); %-0.21;
n_lieblein_hub = lieblein_ncoeff(0,HUBfp.sigma_IGV); %-0.19;

% Optimal incidence angle
i0_tip = Ki_thick_tip * Ki_shape * i010_tip;
i0_mid = Ki_thick_mid * Ki_shape * i010_mid;
i0_hub = Ki_thick_hub * Ki_shape * i010_hub;

TIP.i_optIGV = i0_tip + n_lieblein_tip * TIPfp.theta_IGV;
MID.i_optIGV = i0_mid + n_lieblein_mid * MIDfp.theta_IGV;
HUB.i_optIGV = i0_hub + n_lieblein_hub * HUBfp.theta_IGV;

HUB.gamma_IGV = 0 + HUB.i_optIGV + HUBfp.theta_IGV / 2;
MID.gamma_IGV = 0 + MID.i_optIGV + MIDfp.theta_IGV / 2;
TIP.gamma_IGV = 0 + TIP.i_optIGV + TIPfp.theta_IGV / 2;
