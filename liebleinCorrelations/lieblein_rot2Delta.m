% Outlet deviation angle
% delta = delta0 + ( m * theta ) / sigma^b
% delta0 = delta010 * Kdelta_thick * Kdelta_shape

%for NACA 65
Kdelta_shape_tip = lieblein_deltacorrectionShape(profile{6});%1;
Kdelta_shape_mid = lieblein_deltacorrectionShape(profile{5});%1;
Kdelta_shape_hub = lieblein_deltacorrectionShape(profile{4});%1;

%correction for the thickness
% do graph --> [Kdelta_thick] = liblein_deltathick(t,c)
Kdelta_thick_tip = lieblein_deltacorrectionThick(thick2);%0.75;
Kdelta_thick_mid = lieblein_deltacorrectionThick(thick2);%0.75;
Kdelta_thick_hub = lieblein_deltacorrectionThick(thick2);%0.75;

%Deviation angle with reference thickness
% do graph --> [delta] = lieblein_delta010(beta1,sigma)
delta010_tip = lieblein_delta010(TIPfp.beta3,TIP.sigma2); %1.6;
delta010_mid = lieblein_delta010(MIDfp.beta3,MID.sigma2); %1.5;
delta010_hub = lieblein_delta010(HUBfp.beta3,HUB.sigma2); %1.3;

%coefficient m
% do graph --> [m] = lieblein_mcoeff(beta1,'profile')
m_lieblein_tip = lieblein_mcoeff(TIPfp.beta3,profile{6});%0.26;
m_lieblein_mid = lieblein_mcoeff(MIDfp.beta3,profile{5});%0.25;
m_lieblein_hub = lieblein_mcoeff(HUBfp.beta3,profile{4});%0.23;

%exponent b
% do graph --> [b] = lieblein_bcoeff(beta1)
b_lieblein_tip = lieblein_bcoeff(TIPfp.beta3);%0.68;
b_lieblein_mid = lieblein_bcoeff(MIDfp.beta3);%0.76;
b_lieblein_hub = lieblein_bcoeff(HUBfp.beta3);%0.81;

% optimal deviation angle
delta0_tip = Kdelta_thick_tip * Kdelta_shape_tip * delta010_tip;
delta0_mid = Kdelta_thick_mid * Kdelta_shape_mid * delta010_mid;
delta0_hub = Kdelta_thick_hub * Kdelta_shape_hub * delta010_hub;

TIPfp.delta2 = delta0_tip + m_lieblein_tip * TIP.theta2 / (TIP.sigma2^b_lieblein_tip);
MIDfp.delta2 = delta0_mid + m_lieblein_mid * MID.theta2 / (MID.sigma2^b_lieblein_mid);
HUBfp.delta2 = delta0_hub + m_lieblein_hub * HUB.theta2 / (HUB.sigma2^b_lieblein_hub);
