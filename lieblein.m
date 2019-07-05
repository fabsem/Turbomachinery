%% Geometrical Angle
% i = i0 + n*theta
% i0 = Kth * Kshape * i0.10

% DO GRAPHS --> [i010] = liblein_incidence010(beta1,sigma)
i010_tip = lieblein_incidence010(beta1_tip,sigma_tip);%4;
i010_mid = lieblein_incidence010(beta1_mid,sigma_mid);%4;
i010_hub = lieblein_incidence010(beta1_hub,sigma_hub);%4.2;

%correction coefficient for the shape
Ki_shape = lieblein_icorrectionShape('naca65'); %1; %for NACA 65

%correction for the thickness
% do graph --> [K_thick] = liblein_ithick(t,c)
%[Ki_thick] = lieblein_icorrectionThick(t_c)
Ki_thick_tip = lieblein_icorrectionThick(th); %0.9;
Ki_thick_mid = lieblein_icorrectionThick(th); %0.9;
Ki_thick_hub = lieblein_icorrectionThick(th); %0.9;

% n coefficient
% do graph --> [n_lieblein] = lieblein_ncoeff(beta1,sigma)
n_lieblein_tip = lieblein_ncoeff(beta1_tip,sigma_tip); %-0.26;
n_lieblein_mid = lieblein_ncoeff(beta1_mid,sigma_mid); %-0.21;
n_lieblein_hub = lieblein_ncoeff(beta1_hub,sigma_hub); %-0.19;

% Optimal incidence angle
i0_tip = Ki_thick_tip * Ki_shape * i010_tip;
i0_mid = Ki_thick_mid * Ki_shape * i010_mid;
i0_hub = Ki_thick_hub * Ki_shape * i010_hub;

i_opt_tip = i0_tip + n_lieblein_tip * theta_tip;
i_opt_mid = i0_mid + n_lieblein_mid * theta_mid;
i_opt_hub = i0_hub + n_lieblein_hub * theta_hub;

% Outlet deviation angle
% delta = delta0 + ( m * theta ) / sigma^b
% delta0 = delta010 * Kdelta_thick * Kdelta_shape

%for NACA 65
Kdelta_shape_tip = lieblein_deltacorrectionShape('naca65');%1;
Kdelta_shape_mid = lieblein_deltacorrectionShape('naca65');%1;
Kdelta_shape_hub = lieblein_deltacorrectionShape('naca65');%1;

%correction for the thickness
% do graph --> [Kdelta_thick] = liblein_deltathick(t,c)
Kdelta_thick_tip = lieblein_deltacorrectionThick(th);%0.75;
Kdelta_thick_mid = lieblein_deltacorrectionThick(th);%0.75;
Kdelta_thick_hub = lieblein_deltacorrectionThick(th);%0.75;

%Deviation angle with reference thickness
% do graph --> [delta] = lieblein_delta010(beta1,sigma)
delta010_tip = lieblein_delta010(beta1_tip,sigma_tip); %1.6;
delta010_mid = lieblein_delta010(beta1_mid,sigma_mid); %1.5;
delta010_hub = lieblein_delta010(beta1_hub,sigma_hub); %1.3;

%coefficient m
% do graph --> [m] = lieblein_mcoeff(beta1,'profile')
m_lieblein_tip = lieblein_mcoeff(beta1_tip,'naca65');%0.26;
m_lieblein_mid = lieblein_mcoeff(beta1_mid,'naca65');%0.25;
m_lieblein_hub = lieblein_mcoeff(beta1_hub,'naca65');%0.23;

%exponent b
% do graph --> [b] = lieblein_bcoeff(beta1)
b_lieblein_tip = lieblein_bcoeff(beta1_tip);%0.68;
b_lieblein_mid = lieblein_bcoeff(beta1_mid);%0.76;
b_lieblein_hub = lieblein_bcoeff(beta1_hub);%0.81;

% optimal deviation angle
delta0_tip = Kdelta_thick_tip * Kdelta_shape_tip * delta010_tip;
delta0_mid = Kdelta_thick_mid * Kdelta_shape_mid * delta010_mid;
delta0_hub = Kdelta_thick_hub * Kdelta_shape_hub * delta010_hub;

delta_tip = delta0_tip + m_lieblein_tip * theta_tip / (sigma_tip^b_lieblein_tip);
delta_mid = delta0_mid + m_lieblein_mid * theta_mid / (sigma_mid^b_lieblein_mid);
delta_hub = delta0_hub + m_lieblein_hub * theta_hub / (sigma_hub^b_lieblein_hub);


%%  Flow deflection according to Lieblein
%epsilon_tip = theta_tip + i_opt_tip - delta_tip;
%epsilon_mid = theta_mid + i_opt_mid - delta_mid;
%epsilon_hub = theta_hub + i_opt_hub - delta_hub;

theta_new_hub = epsilon_hub - i_opt_hub + delta_hub;
theta_new_mid = epsilon_mid - i_opt_mid + delta_mid;
theta_new_tip = epsilon_tip - i_opt_tip + delta_tip;


    err(1) = abs((theta_new_hub - theta_hub)/theta_hub);
    err(2) = abs((theta_new_mid - theta_mid)/theta_mid);
    err(3) = abs((theta_new_tip - theta_tip)/theta_tip);

    err = max(err);

    theta_hub = theta_new_hub;
    theta_mid = theta_new_mid;
    theta_tip = theta_new_tip;

    i = i + 1;
