%% Geometrical Angle
% i = i0 + n*theta
% i0 = Kth * Kshape * i0.10

profile_IGV = 'naca65';
thick_IGV = 0.1;

% DO GRAPHS --> [i010] = liblein_incidence010(beta1,sigma)
i010_tip = lieblein_incidence010(0,TIPfp.sigma_IGV);%4;
i010_mid = lieblein_incidence010(0,MIDfp.sigma_IGV);%4;
i010_hub = lieblein_incidence010(0,HUBfp.sigma_IGV);%4.2;

%correction coefficient for the shape
Ki_shape = lieblein_icorrectionShape(profile_IGV); %1; %for NACA 65

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

TIPfp.i_IGV = i0_tip + n_lieblein_tip * TIPfp.theta_IGV;
MIDfp.i_IGV = i0_mid + n_lieblein_mid * MIDfp.theta_IGV;
HUBfp.i_IGV = i0_hub + n_lieblein_hub * HUBfp.theta_IGV;

% Outlet deviation angle
% delta = delta0 + ( m * theta ) / sigma^b
% delta0 = delta010 * Kdelta_thick * Kdelta_shape

%for NACA 65
Kdelta_shape_tip = lieblein_deltacorrectionShape(profile_IGV);%1;
Kdelta_shape_mid = lieblein_deltacorrectionShape(profile_IGV);%1;
Kdelta_shape_hub = lieblein_deltacorrectionShape(profile_IGV);%1;

%correction for the thickness
% do graph --> [Kdelta_thick] = liblein_deltathick(t,c)
Kdelta_thick_tip = lieblein_deltacorrectionThick(thick_IGV);%0.75;
Kdelta_thick_mid = lieblein_deltacorrectionThick(thick_IGV);%0.75;
Kdelta_thick_hub = lieblein_deltacorrectionThick(thick_IGV);%0.75;

%Deviation angle with reference thickness
% do graph --> [delta] = lieblein_delta010(beta1,sigma)
delta010_tip = lieblein_delta010(0,TIPfp.sigma_IGV); %1.6;
delta010_mid = lieblein_delta010(0,MIDfp.sigma_IGV); %1.5;
delta010_hub = lieblein_delta010(0,HUBfp.sigma_IGV); %1.3;

%coefficient m
% do graph --> [m] = lieblein_mcoeff(beta1,'profile')
m_lieblein_tip = lieblein_mcoeff(0,profile_IGV);%0.26;
m_lieblein_mid = lieblein_mcoeff(0,profile_IGV);%0.25;
m_lieblein_hub = lieblein_mcoeff(0,profile_IGV);%0.23;

%exponent b
% do graph --> [b] = lieblein_bcoeff(beta1)
b_lieblein_tip = lieblein_bcoeff(0);%0.68;
b_lieblein_mid = lieblein_bcoeff(0);%0.76;
b_lieblein_hub = lieblein_bcoeff(0);%0.81;

% optimal deviation angle
delta0_tip = Kdelta_thick_tip * Kdelta_shape_tip * delta010_tip;
delta0_mid = Kdelta_thick_mid * Kdelta_shape_mid * delta010_mid;
delta0_hub = Kdelta_thick_hub * Kdelta_shape_hub * delta010_hub;

TIPfp.delta_IGV = delta0_tip + m_lieblein_tip * TIPfp.theta_IGV / (TIPfp.sigma_IGV^b_lieblein_tip);
MIDfp.delta_IGV = delta0_mid + m_lieblein_mid * MIDfp.theta_IGV / (MIDfp.sigma_IGV^b_lieblein_mid);
HUBfp.delta_IGV = delta0_hub + m_lieblein_hub * HUBfp.theta_IGV / (HUBfp.sigma_IGV^b_lieblein_hub);


%%  Flow deflection according to Lieblein
%epsilon_tip = theta_tip + i_opt_tip - delta_tip;
%epsilon_mid = theta_mid + i_opt_mid - delta_mid;
%epsilon_hub = theta_hub + i_opt_hub - delta_hub;

theta_new_hub = HUBfp.epsilon_IGV - HUBfp.i_IGV + HUBfp.delta_IGV;
theta_new_mid = MIDfp.epsilon_IGV - MIDfp.i_IGV + MIDfp.delta_IGV;
theta_new_tip = TIPfp.epsilon_IGV - TIPfp.i_IGV+ TIPfp.delta_IGV;


    err(1) = abs((theta_new_hub - HUBfp.theta_IGV)/HUBfp.theta_IGV);
    err(2) = abs((theta_new_mid - MIDfp.theta_IGV)/MIDfp.theta_IGV);
    err(3) = abs((theta_new_tip - TIPfp.theta_IGV)/TIPfp.theta_IGV);

    err = max(err);

    HUBfp.theta_IGV = theta_new_hub;
    MIDfp.theta_IGV = theta_new_mid;
    TIPfp.theta_IGV = theta_new_tip;

    i = i + 1;
