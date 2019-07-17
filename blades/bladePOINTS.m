% Run inside main script

%print with stagger
stagger = 1;

%Rotor1
[x1,y1] = bladeMaker('DCA',HUB.theta1,HUB.c1,stagger,HUB.gamma1,'Hub - Rotor 1',1);
[x2,y2] = bladeMaker('DCA',MID.theta1,MID.c1,stagger,MID.gamma1,'Mid - Rotor 1',1);
[x3,y3] = bladeMaker('DCA',TIP.theta1,TIP.c1,stagger,TIP.gamma1,'Tip - Rotor 1',1);

writeGeomFile(y1,[],x1,'hub_R1')
writeGeomFile(y2,[],x2,'mid_R1')
writeGeomFile(y3,[],x3,'tip_R1')

%Rotor2
[x4,y4] = bladeMaker('DCA',HUB.theta2,HUB.c2,stagger,HUB.gamma2,'Hub - Rotor 2',2);
[x5,y5] = bladeMaker('DCA',MID.theta2,MID.c2,stagger,MID.gamma2,'Mid - Rotor 2',2);
[x6,y6] = bladeMaker('DCA',TIP.theta2,TIP.c2,stagger,TIP.gamma2,'Tip - Rotor 2',2);

writeGeomFile(y4,[],x4,'hub_R2')
writeGeomFile(y5,[],x5,'mid_R2')
writeGeomFile(y6,[],x6,'tip_R2')

% IGV On Design
[x7,y7] = bladeMaker('naca65',HUBfp.theta_IGV,HUBfp.c_IGV,stagger,HUB.gamma_IGV,'Hub IGV - On Design',[]);
[x8,y8] = bladeMaker('naca65',MIDfp.theta_IGV,MIDfp.c_IGV,stagger,MID.gamma_IGV,'Mid IGV - On Design',[]);
[x9,y9] = bladeMaker('naca65',TIPfp.theta_IGV,TIPfp.c_IGV,stagger,TIP.gamma_IGV,'Tip IGV - On Design',[]);

writeGeomFile(y7,[],x7,'hub_IGV_on')
writeGeomFile(y8,[],x8,'mid_IGV_on')
writeGeomFile(y9,[],x9,'tip_IGV_on')

% IGV Off Design
[x10,y10] = bladeMaker('naca65',HUBfp.theta_IGV,HUBfp.c_IGV,stagger,HUBfp.gamma_IGV,'Hub IGV - Off Design',[]);
[x11,y11] = bladeMaker('naca65',MIDfp.theta_IGV,MIDfp.c_IGV,stagger,MIDfp.gamma_IGV,'Mid IGV - Off Design',[]);
[x12,y12] = bladeMaker('naca65',TIPfp.theta_IGV,TIPfp.c_IGV,stagger,TIPfp.gamma_IGV,'Tip IGV - Off Design',[]);

writeGeomFile(y10,[],x10,'hub_IGV_off')
writeGeomFile(y11,[],x11,'mid_IGV_off')
writeGeomFile(y12,[],x12,'tip_IGV_off')