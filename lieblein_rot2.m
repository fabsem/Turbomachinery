%% Geometrical Angle
% i = i0 + n*theta
% i0 = Kth * Kshape * i0.10

profile = 'naca65';

% DO GRAPHS --> [i010] = liblein_incidence010(beta1,sigma)
i010_tip = lieblein_incidence010(TIP.beta3,TIP.sigma2);%4;
i010_mid = lieblein_incidence010(MID.beta3,MID.sigma2);%4;
i010_hub = lieblein_incidence010(HUB.beta3,HUB.sigma2);%4.2;

%correction coefficient for the shape
Ki_shape = lieblein_icorrectionShape(profile); %1; %for NACA 65

%correction for the thickness
% do graph --> [K_thick] = liblein_ithick(t,c)
%[Ki_thick] = lieblein_icorrectionThick(t_c)
Ki_thick_tip = lieblein_icorrectionThick(thick2); %0.9;
Ki_thick_mid = lieblein_icorrectionThick(thick2); %0.9;
Ki_thick_hub = lieblein_icorrectionThick(thick2); %0.9;

% n coefficient
% do graph --> [n_lieblein] = lieblein_ncoeff(beta1,sigma)
n_lieblein_tip = lieblein_ncoeff(TIP.beta3,TIP.sigma2); %-0.26;
n_lieblein_mid = lieblein_ncoeff(MID.beta3,MID.sigma2); %-0.21;
n_lieblein_hub = lieblein_ncoeff(HUB.beta3,HUB.sigma2); %-0.19;

% Optimal incidence angle
i0_tip = Ki_thick_tip * Ki_shape * i010_tip;
i0_mid = Ki_thick_mid * Ki_shape * i010_mid;
i0_hub = Ki_thick_hub * Ki_shape * i010_hub;

i_opt_tip = i0_tip + n_lieblein_tip * TIP.theta2;
i_opt_mid = i0_mid + n_lieblein_mid * MID.theta2;
i_opt_hub = i0_hub + n_lieblein_hub * HUB.theta2;

% Outlet deviation angle
% delta = delta0 + ( m * theta ) / sigma^b
% delta0 = delta010 * Kdelta_thick * Kdelta_shape

%for NACA 65
Kdelta_shape_tip = lieblein_deltacorrectionShape(profile);%1;
Kdelta_shape_mid = lieblein_deltacorrectionShape(profile);%1;
Kdelta_shape_hub = lieblein_deltacorrectionShape(profile);%1;

%correction for the thickness
% do graph --> [Kdelta_thick] = liblein_deltathick(t,c)
Kdelta_thick_tip = lieblein_deltacorrectionThick(thick2);%0.75;
Kdelta_thick_mid = lieblein_deltacorrectionThick(thick2);%0.75;
Kdelta_thick_hub = lieblein_deltacorrectionThick(thick2);%0.75;

%Deviation angle with reference thickness
% do graph --> [delta] = lieblein_delta010(beta1,sigma)
delta010_tip = lieblein_delta010(TIP.beta3,TIP.sigma2); %1.6;
delta010_mid = lieblein_delta010(MID.beta3,MID.sigma2); %1.5;
delta010_hub = lieblein_delta010(HUB.beta3,HUB.sigma2); %1.3;

%coefficient m
% do graph --> [m] = lieblein_mcoeff(beta1,'profile')
m_lieblein_tip = lieblein_mcoeff(TIP.beta3,profile);%0.26;
m_lieblein_mid = lieblein_mcoeff(MID.beta3,profile);%0.25;
m_lieblein_hub = lieblein_mcoeff(HUB.beta3,profile);%0.23;

%exponent b
% do graph --> [b] = lieblein_bcoeff(beta1)
b_lieblein_tip = lieblein_bcoeff(TIP.beta3);%0.68;
b_lieblein_mid = lieblein_bcoeff(MID.beta3);%0.76;
b_lieblein_hub = lieblein_bcoeff(HUB.beta3);%0.81;

% optimal deviation angle
delta0_tip = Kdelta_thick_tip * Kdelta_shape_tip * delta010_tip;
delta0_mid = Kdelta_thick_mid * Kdelta_shape_mid * delta010_mid;
delta0_hub = Kdelta_thick_hub * Kdelta_shape_hub * delta010_hub;

TIP.delta2 = delta0_tip + m_lieblein_tip * TIP.theta2 / (TIP.sigma2^b_lieblein_tip);
MID.delta2 = delta0_mid + m_lieblein_mid * MID.theta2 / (MID.sigma2^b_lieblein_mid);
HUB.delta2 = delta0_hub + m_lieblein_hub * HUB.theta2 / (HUB.sigma2^b_lieblein_hub);


%%  Flow deflection according to Lieblein
%epsilon_tip = theta_tip + i_opt_tip - delta_tip;
%epsilon_mid = theta_mid + i_opt_mid - delta_mid;
%epsilon_hub = theta_hub + i_opt_hub - delta_hub;

theta_new_hub = HUB.epsilon2 - i_opt_hub + HUB.delta2;
theta_new_mid = MID.epsilon2 - i_opt_mid + MID.delta2;
theta_new_tip = TIP.epsilon2 - i_opt_tip + TIP.delta2;


    err(1) = abs((theta_new_hub - HUB.theta2)/HUB.theta2);
    err(2) = abs((theta_new_mid - MID.theta2)/MID.theta2);
    err(3) = abs((theta_new_tip - TIP.theta2)/TIP.theta2);

    err = max(err);

    HUB.theta2 = theta_new_hub;
    MID.theta2 = theta_new_mid;
    TIP.theta2 = theta_new_tip;

    i = i + 1;
