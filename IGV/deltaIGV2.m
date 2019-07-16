function [hub_delta2 mid_delta2 tip_delta2] = deltaIGV2 (HUBfp,MIDfp,TIPfp,HUB,MID,TIP,profile,thick2)
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
delta010_tip = lieblein_delta010(TIPfp.beta2,TIP.sigma2); %1.6;
delta010_mid = lieblein_delta010(MIDfp.beta2,MID.sigma2); %1.5;
delta010_hub = lieblein_delta010(HUBfp.beta2,HUB.sigma2); %1.3;

%coefficient m
% do graph --> [m] = lieblein_mcoeff(beta1,'profile')
m_lieblein_tip = lieblein_mcoeff(TIPfp.beta2,profile);%0.26;
m_lieblein_mid = lieblein_mcoeff(MIDfp.beta2,profile);%0.25;
m_lieblein_hub = lieblein_mcoeff(HUBfp.beta2,profile);%0.23;

%exponent b
% do graph --> [b] = lieblein_bcoeff(beta1)
b_lieblein_tip = lieblein_bcoeff(TIPfp.beta2);%0.68;
b_lieblein_mid = lieblein_bcoeff(MIDfp.beta2);%0.76;
b_lieblein_hub = lieblein_bcoeff(HUBfp.beta2);%0.81;

% optimal deviation angle
delta0_tip = Kdelta_thick_tip * Kdelta_shape_tip * delta010_tip;
delta0_mid = Kdelta_thick_mid * Kdelta_shape_mid * delta010_mid;
delta0_hub = Kdelta_thick_hub * Kdelta_shape_hub * delta010_hub;

tip_delta2 = delta0_tip + m_lieblein_tip * TIP.theta2 / (TIP.sigma1^b_lieblein_tip);
mid_delta2 = delta0_mid + m_lieblein_mid * MID.theta2 / (MID.sigma1^b_lieblein_mid);
hub_delta2 = delta0_hub + m_lieblein_hub * HUB.theta2 / (HUB.sigma1^b_lieblein_hub);
