clear all
close all
clc

c = 1;

theta = 30;
tb = 0.08 * c;
r0 = 0.01 * c;
theta_span = -theta:theta;

Rc = c/2 / sind(theta/2);
yc = -Rc * cosd(theta/2);
y0 = c/2 * tand(theta/4);
x = linspace(-0.5,0.5);
y = yc + sqrt(Rc^2 - x.^2);


%Cl0 = tand(theta/4) / 0.1103;
%Cl = y0 / c * Cl0

% d
du = y0 + tb/2 - r0 * sind(theta/2);

%Upper surface
Ru = (du^2 - r0^2 + (c/2 - r0 * cosd(theta/2))^2) / (2 * (du - r0));
ycu = y0 + tb/2 - Ru;

xu = x;
yu = ycu + sqrt(Ru^2 - x.^2);

%Lower surface
dl = y0 - tb/2 + r0 * sind(theta/2);
Rl = (dl^2 - r0^2 + (c/2 + r0 * cosd(theta/2))^2) / (2 * (dl + r0));
ycl = y0 - tb/2 - Rl;

xl = x;
yl = ycl + sqrt(Rl^2 - x.^2);

%smooth edges
y_edge0_dx = r0 * sind(theta/2);
x_edge0_dx = + (c/2 - r0 * cosd(theta/2));

theta_r0 = 0:360;%-100:65;
x_edge_dx = r0 * cosd(theta_r0) + x_edge0_dx;
y_edge_dx = r0 * sind(theta_r0) + y_edge0_dx;

y_edge0_sx = r0 * sind(theta/2);
x_edge0_sx = - (c/2 - r0 * cosd(theta/2));

theta_r0 = 0:360;%115:260;
x_edge_sx = r0 * cosd(theta_r0) + x_edge0_sx;
y_edge_sx = r0 * sind(theta_r0) + y_edge0_sx;

plot(x,y,'--')
hold on
plot(xu,yu,'linewidth',3)
plot(xl,yl,'linewidth',3)
plot(x_edge_dx,y_edge_dx)
plot(x_edge_sx,y_edge_sx)
title(['DCA Profile - c = ', num2str(c),', t = ',num2str(tb),'c'])
axis equal

%% Output vector
x = [flip(xu), xl]';
y = [flip(yu), yl]';