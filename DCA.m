clear all
close all

c = 1;

theta = 30;
tb = 0.1 * c;
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

plot(x,y)
hold on
plot(xu,yu)
plot(xl,yl)
axis equal