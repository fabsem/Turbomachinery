function printBlades(profile,gamma_hub,gamma_mid,gamma_tip,c_hub,c_mid,c_tip,b,theta_hub,theta_mid,theta_tip)



% select shape
if strcmp(profile,'naca0012')
  gamma_hub = -26;
  gamma_mid = -41;
  gamma_tip = -50.7;

  c = 0.13;

  b = 0.1736;
x = [1;0.999416000000000;0.997666000000000;0.994753000000000;0.990685000000000;0.985471000000000;0.979123000000000;0.971656000000000;0.963087000000000;0.953437000000000;0.942728000000000;0.930985000000000;0.918235000000000;0.904509000000000;0.889837000000000;0.874255000000000;0.857800000000000;0.840508000000000;0.822421000000000;0.803581000000000;0.784032000000000;0.763820000000000;0.742992000000000;0.721596000000000;0.699682000000000;0.677303000000000;0.654509000000000;0.631354000000000;0.607892000000000;0.584179000000000;0.560268000000000;0.536217000000000;0.512082000000000;0.487918000000000;0.463783000000000;0.439732000000000;0.415822000000000;0.392108000000000;0.368646000000000;0.345492000000000;0.322698000000000;0.300318000000000;0.278404000000000;0.257008000000000;0.236180000000000;0.215968000000000;0.196419000000000;0.177579000000000;0.159492000000000;0.142201000000000;0.125745000000000;0.110163000000000;0.0954920000000000;0.0817650000000000;0.0690150000000000;0.0572720000000000;0.0465630000000000;0.0369130000000000;0.0283440000000000;0.0208770000000000;0.0145290000000000;0.00931500000000000;0.00524700000000000;0.00233400000000000;0.000584000000000000;0;0.000584000000000000;0.00233400000000000;0.00524700000000000;0.00931500000000000;0.0145290000000000;0.0208770000000000;0.0283440000000000;0.0369130000000000;0.0465630000000000;0.0572720000000000;0.0690150000000000;0.0817650000000000;0.0954920000000000;0.110163000000000;0.125745000000000;0.142201000000000;0.159492000000000;0.177579000000000;0.196419000000000;0.215968000000000;0.236180000000000;0.257008000000000;0.278404000000000;0.300318000000000;0.322698000000000;0.345492000000000;0.368646000000000;0.392108000000000;0.415822000000000;0.439732000000000;0.463783000000000;0.487918000000000;0.512082000000000;0.536217000000000;0.560268000000000;0.584179000000000;0.607892000000000;0.631354000000000;0.654509000000000;0.677303000000000;0.699682000000000;0.721596000000000;0.742992000000000;0.763820000000000;0.784032000000000;0.803581000000000;0.822421000000000;0.840508000000000;0.857800000000000;0.874255000000000;0.889837000000000;0.904509000000000;0.918235000000000;0.930985000000000;0.942728000000000;0.953437000000000;0.963087000000000;0.971656000000000;0.979123000000000;0.985471000000000;0.990685000000000;0.994753000000000;0.997666000000000;0.999416000000000;1];
y = [0.00126000000000000;0.00134200000000000;0.00158700000000000;0.00199400000000000;0.00256000000000000;0.00328000000000000;0.00415200000000000;0.00516900000000000;0.00632400000000000;0.00761100000000000;0.00902200000000000;0.0105490000000000;0.0121820000000000;0.0139140000000000;0.0157350000000000;0.0176350000000000;0.0196050000000000;0.0216350000000000;0.0237140000000000;0.0258340000000000;0.0279830000000000;0.0301520000000000;0.0323290000000000;0.0345060000000000;0.0366700000000000;0.0388110000000000;0.0409170000000000;0.0429780000000000;0.0449800000000000;0.0469120000000000;0.0487620000000000;0.0505160000000000;0.0521620000000000;0.0536870000000000;0.0550770000000000;0.0563200000000000;0.0574030000000000;0.0583140000000000;0.0590420000000000;0.0595750000000000;0.0599030000000000;0.0600170000000000;0.0599100000000000;0.0595760000000000;0.0590080000000000;0.0582050000000000;0.0571640000000000;0.0558860000000000;0.0543720000000000;0.0526250000000000;0.0506510000000000;0.0484570000000000;0.0460490000000000;0.0434370000000000;0.0406310000000000;0.0376410000000000;0.0344790000000000;0.0311560000000000;0.0276830000000000;0.0240710000000000;0.0203300000000000;0.0164710000000000;0.0125010000000000;0.00842900000000000;0.00426000000000000;0;-0.00426000000000000;-0.00842900000000000;-0.0125010000000000;-0.0164710000000000;-0.0203300000000000;-0.0240710000000000;-0.0276830000000000;-0.0311560000000000;-0.0344790000000000;-0.0376410000000000;-0.0406310000000000;-0.0434370000000000;-0.0460490000000000;-0.0484570000000000;-0.0506510000000000;-0.0526250000000000;-0.0543720000000000;-0.0558860000000000;-0.0571640000000000;-0.0582050000000000;-0.0590080000000000;-0.0595760000000000;-0.0599100000000000;-0.0600170000000000;-0.0599030000000000;-0.0595750000000000;-0.0590420000000000;-0.0583140000000000;-0.0574030000000000;-0.0563200000000000;-0.0550770000000000;-0.0536870000000000;-0.0521620000000000;-0.0505160000000000;-0.0487620000000000;-0.0469120000000000;-0.0449800000000000;-0.0429780000000000;-0.0409170000000000;-0.0388110000000000;-0.0366700000000000;-0.0345060000000000;-0.0323290000000000;-0.0301520000000000;-0.0279830000000000;-0.0258340000000000;-0.0237140000000000;-0.0216350000000000;-0.0196050000000000;-0.0176350000000000;-0.0157350000000000;-0.0139140000000000;-0.0121820000000000;-0.0105490000000000;-0.00902200000000000;-0.00761100000000000;-0.00632400000000000;-0.00516900000000000;-0.00415200000000000;-0.00328000000000000;-0.00256000000000000;-0.00199400000000000;-0.00158700000000000;-0.00134200000000000;-0.00126000000000000];

elseif strcmp(profile,'DCA')

    % Due to different reference, invert the stagger
    gamma_hub = gamma_hub * -1;
    gamma_mid = gamma_mid * -1;
    gamma_tip = gamma_tip * -1;
    
[x1,y1] = bladeMaker(profile,theta_hub,c_hub);
[x2,y2] = bladeMaker(profile,theta_mid,c_mid);
[x3,y3] = bladeMaker(profile,theta_tip,c_tip);

scale = 1;%/c_hub;

figure
%plot3(linspace(-c_hub/2,c_hub/2) * cosd(gamma_hub), linspace(-c_hub/2,c_hub/2) * sind(gamma_hub), zeros(1,100))
hold on
%plot3(linspace(-c_mid/2,c_mid/2) * cosd(gamma_mid), linspace(-c_mid/2,c_mid/2) * sind(gamma_mid), ones(1,100) * b/2)
%plot3(linspace(-c_tip/2,c_tip/2) * cosd(gamma_tip), linspace(-c_tip/2,c_tip/2) * sind(gamma_tip), ones(1,100) * b)
plot3(x1./scale .* cosd(gamma_hub) - y1./scale .* sind(gamma_hub) ,x1./scale .* sind(gamma_hub) + y1./scale .* cosd(gamma_hub),zeros(1,length(x1)))
plot3(x2./scale .* cosd(gamma_mid) - y2./scale .* sind(gamma_mid) ,x2./scale .* sind(gamma_mid) + y2./scale .* cosd(gamma_mid),ones(1,length(x1)) * b/2)
plot3(x3./scale .* cosd(gamma_tip) - y3./scale .* sind(gamma_tip) ,x3./scale .* sind(gamma_tip) + y3./scale .* cosd(gamma_tip),ones(1,length(x1)) * b)
surface([x1./scale .* cosd(gamma_hub) - y1./scale .* sind(gamma_hub), x2./scale .* cosd(gamma_mid) - y2./scale .* sind(gamma_mid)], [x1./scale .* sind(gamma_hub) + y1./scale .* cosd(gamma_hub), x2./scale .* sind(gamma_mid) + y2./scale .* cosd(gamma_mid)],[zeros(1,length(x1))', ones(1,length(x1))' * b/2])
surface([x2./scale .* cosd(gamma_mid) - y2/scale .* sind(gamma_mid), x3./scale .* cosd(gamma_tip) - y3./scale .* sind(gamma_tip)],[x2./scale .* sind(gamma_mid) + y2./scale .* cosd(gamma_mid),x3./scale .* sind(gamma_tip) + y3./scale .* cosd(gamma_tip)],[ones(1,length(x1))' * b/2,  ones(1,length(x1))' * b])
axis equal
end
end
