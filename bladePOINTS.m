% Run inside main script

%Rotor1
[x1,y1] = bladeMaker('DCA',HUB.theta1,HUB.c1);
[x2,y2] = bladeMaker('DCA',MID.theta1,MID.c1);
[x3,y3] = bladeMaker('DCA',TIP.theta1,TIP.c1);

writeGeomFile(x1,y1,[],'hub_R1')
writeGeomFile(x2,y2,[],'mid_R1')
writeGeomFile(x3,y3,[],'tip_R1')

%Rotor2
[x4,y4] = bladeMaker('DCA',HUB.theta2,HUB.c2);
[x5,y5] = bladeMaker('DCA',MID.theta2,MID.c2);
[x6,y6] = bladeMaker('DCA',TIP.theta2,TIP.c2);

writeGeomFile(x4,y4,[],'hub_R2')
writeGeomFile(x5,y5,[],'mid_R2')
writeGeomFile(x6,y6,[],'tip_R2')
