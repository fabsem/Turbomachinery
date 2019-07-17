% Run inside main script
clear all
close all
clc
load('results_contra');

%Rotor1
[x1,y1] = bladeMaker('DCA',HUB.theta1,HUB.c1);
[x2,y2] = bladeMaker('DCA',MID.theta1,MID.c1);
[x3,y3] = bladeMaker('DCA',TIP.theta1,TIP.c1);

%Add stagger
gamma_hub1 = HUB.gamma1 * -1;
gamma_mid1 = MID.gamma1 * -1;
gamma_tip1 = TIP.gamma1 * -1;

R = @(gamma) [cosd(gamma) sind(gamma);
                -sind(gamma) cosd(gamma)];
 
            
coord_hub1 = (R(gamma_hub1) * [x1,y1]')';
coord_mid1 = (R(gamma_mid1) * [x2,y2]')';
coord_tip1 = (R(gamma_tip1) * [x3,y3]')';

x1 = coord_hub1(:,1); y1 = coord_hub1(:,2);
x2 = coord_mid1(:,1); y2 = coord_mid1(:,2);
x3 = coord_tip1(:,1); y3 = coord_tip1(:,2);

writeGeomFile(x1,y1,[],'hub_R1')
writeGeomFile(x2,y2,[],'mid_R1')
writeGeomFile(x3,y3,[],'tip_R1')

%Rotor2
[x4,y4] = bladeMaker('DCA',HUB.theta2,HUB.c2);
[x5,y5] = bladeMaker('DCA',MID.theta2,MID.c2);
[x6,y6] = bladeMaker('DCA',TIP.theta2,TIP.c2);

%Add stagger
gamma_hub2 = HUB.gamma2 * -1;
gamma_mid2 = MID.gamma2 * -1;
gamma_tip2 = TIP.gamma2 * -1;

R = @(gamma) [cosd(gamma) sind(gamma);
                -sind(gamma) cosd(gamma)];
 
            
coord_hub2 = (R(gamma_hub2) * [x4,y4]')';
coord_mid2 = (R(gamma_mid2) * [x5,y5]')';
coord_tip2 = (R(gamma_tip2) * [x6,y6]')';

x4 = coord_hub2(:,1); y4 = coord_hub2(:,2);
x5 = coord_mid2(:,1); y5 = coord_mid2(:,2);
x6 = coord_tip2(:,1); y6 = coord_tip2(:,2);

writeGeomFile(x4,y4,[],'hub_R2')
writeGeomFile(x5,y5,[],'mid_R2')
writeGeomFile(x6,y6,[],'tip_R2')