%% Geometrical Angle
% i = i0 + n*theta
% i0 = Kth * Kshape * i0.10

profile = 'DCA';

% DO GRAPHS --> [i010] = liblein_incidence010(beta1,sigma)
i010_tip = lieblein_incidence010(TIP.beta1,TIP.sigma1);%4;
i010_mid = lieblein_incidence010(MID.beta1,MID.sigma1);%4;
i010_hub = lieblein_incidence010(HUB.beta1,HUB.sigma1);%4.2;

%correction coefficient for the shape
Ki_shape = lieblein_icorrectionShape(profile); %1; %for NACA 65

%correction for the thickness
% do graph --> [K_thick] = liblein_ithick(t,c)
%[Ki_thick] = lieblein_icorrectionThick(t_c)
Ki_thick_tip = lieblein_icorrectionThick(thick1); %0.9;
Ki_thick_mid = lieblein_icorrectionThick(thick1); %0.9;
Ki_thick_hub = lieblein_icorrectionThick(thick1); %0.9;

% n coefficient
% do graph --> [n_lieblein] = lieblein_ncoeff(beta1,sigma)
n_lieblein_tip = lieblein_ncoeff(TIP.beta1,TIP.sigma1); %-0.26;
n_lieblein_mid = lieblein_ncoeff(MID.beta1,MID.sigma1); %-0.21;
n_lieblein_hub = lieblein_ncoeff(HUB.beta1,HUB.sigma1); %-0.19;

% Optimal incidence angle
i0_tip = Ki_thick_tip * Ki_shape * i010_tip;
i0_mid = Ki_thick_mid * Ki_shape * i010_mid;
i0_hub = Ki_thick_hub * Ki_shape * i010_hub;

TIP.i_opt1 = i0_tip + n_lieblein_tip * TIP.theta1;
MID.i_opt1 = i0_mid + n_lieblein_mid * MID.theta1;
HUB.i_opt1 = i0_hub + n_lieblein_hub * HUB.theta1;

% Outlet deviation angle
% delta = delta0 + ( m * theta ) / sigma^b
% delta0 = delta010 * Kdelta_thick * Kdelta_shape

%for NACA 65
Kdelta_shape_tip = lieblein_deltacorrectionShape(profile);%1;
Kdelta_shape_mid = lieblein_deltacorrectionShape(profile);%1;
Kdelta_shape_hub = lieblein_deltacorrectionShape(profile);%1;

%correction for the thickness
% do graph --> [Kdelta_thick] = liblein_deltathick(t,c)
Kdelta_thick_tip = lieblein_deltacorrectionThick(thick1);%0.75;
Kdelta_thick_mid = lieblein_deltacorrectionThick(thick1);%0.75;
Kdelta_thick_hub = lieblein_deltacorrectionThick(thick1);%0.75;

%Deviation angle with reference thickness
% do graph --> [delta] = lieblein_delta010(beta1,sigma)
delta010_tip = lieblein_delta010(TIP.beta1,TIP.sigma1); %1.6;
delta010_mid = lieblein_delta010(MID.beta1,MID.sigma1); %1.5;
delta010_hub = lieblein_delta010(HUB.beta1,HUB.sigma1); %1.3;

%coefficient m
% do graph --> [m] = lieblein_mcoeff(beta1,'profile')
m_lieblein_tip = lieblein_mcoeff(TIP.beta1,profile);%0.26;
m_lieblein_mid = lieblein_mcoeff(MID.beta1,profile);%0.25;
m_lieblein_hub = lieblein_mcoeff(HUB.beta1,profile);%0.23;

%exponent b
% do graph --> [b] = lieblein_bcoeff(beta1)
b_lieblein_tip = lieblein_bcoeff(TIP.beta1);%0.68;
b_lieblein_mid = lieblein_bcoeff(MID.beta1);%0.76;
b_lieblein_hub = lieblein_bcoeff(HUB.beta1);%0.81;

% optimal deviation angle
delta0_tip = Kdelta_thick_tip * Kdelta_shape_tip * delta010_tip;
delta0_mid = Kdelta_thick_mid * Kdelta_shape_mid * delta010_mid;
delta0_hub = Kdelta_thick_hub * Kdelta_shape_hub * delta010_hub;

TIP.delta1 = delta0_tip + m_lieblein_tip * TIP.theta1 / (TIP.sigma1^b_lieblein_tip);
MID.delta1 = delta0_mid + m_lieblein_mid * MID.theta1 / (MID.sigma1^b_lieblein_mid);
HUB.delta1 = delta0_hub + m_lieblein_hub * HUB.theta1 / (HUB.sigma1^b_lieblein_hub);


%%  Flow deflection according to Lieblein
%epsilon_tip = theta_tip + i_opt_tip - delta_tip;
%epsilon_mid = theta_mid + i_opt_mid - delta_mid;
%epsilon_hub = theta_hub + i_opt_hub - delta_hub;

theta_new_hub = HUB.epsilon1 - HUB.i_opt1 + HUB.delta1;
theta_new_mid = MID.epsilon1 - MID.i_opt1 + MID.delta1;
theta_new_tip = TIP.epsilon1 - TIP.i_opt1 + TIP.delta1;


    err(1) = abs((theta_new_hub - HUB.theta1)/HUB.theta1);
    err(2) = abs((theta_new_mid - MID.theta1)/MID.theta1);
    err(3) = abs((theta_new_tip - TIP.theta1)/TIP.theta1);

    err = max(err);

    HUB.theta1 = theta_new_hub;
    MID.theta1 = theta_new_mid;
    TIP.theta1 = theta_new_tip;

    i = i + 1;
