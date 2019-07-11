function [x_profile,y_profile] = bladeMaker(name,theta,c)

theta_blade = abs(theta);

if strcmp(name,'naca65')
%constants
x_c=[0 0.5 0.75 1.25 2.5 5 7.5 10 15 20 25 30 35 40 45 50 ...
     55 60 65 70 75 80 85 90 95 100];% chord percenage x/c
t_c=[0 0.722 0.932 1.169 1.574 2.177 2.647 3.040 3.666 4.143 4.503 4.76 4.924 4.996 4.963 4.812 ...
     4.53 4.146 3.682 3.156 2.584 1.987 1.385 0.81 0.306 0];% half thickness t/c
y_camber=[0 0.25 0.35 0.535 0.93 1.58 2.12 2.585 3.365 3.98 4.475 4.86 5.15 5.355 5.475 5.515 ...
          5.475 5.355 5.15 4.86 4.475 3.98 3.365 2.585 1.58 0];%camber line y/c

cl0=theta_blade/25;
%% equations
x=c*x_c;
t=c*t_c;
y=c*y_camber*cl0;
yu=y+t;
yl=y-t;

figure
plot(x,y,'k--')
hold on
axis equal
plot(x,yu,'k','linewidth',2)
plot(x,yl,'k','linewidth',2)
xlabel('x')
ylabel('y')
title(['NACA 65(',num2str(cl0 * 10),')']);

%output vector
x_profile = [flip(x(end:2)) x]';
y_profile = [flip(y(end:2)) y]';

elseif strcmp(name,'DCA')

  tb = 0.1 * c;
  r0 = 0; %0.01 * c; %POSSIAMO FARE RACCORDO CON SPLINE TRA LE 3 LINEE!!
  theta_span = -theta_blade:theta_blade;

  Rc = c/2 / sind(theta_blade/2);
  yc = -Rc * cosd(theta_blade/2);
  y0 = c/2 * tand(theta_blade/4);
  x = linspace(-0.5 * c,0.5 * c);
  y = yc + sqrt(Rc^2 - x.^2);


  % d
  du = y0 + tb/2 - r0 * sind(theta_blade/2);

  %Upper surface
  Ru = (du^2 - r0^2 + (c/2 - r0 * cosd(theta_blade/2))^2) / (2 * (du - r0));
  ycu = y0 + tb/2 - Ru;

  xu = x;
  yu = ycu + sqrt(Ru^2 - x.^2);

  %Lower surface
  dl = y0 - tb/2 + r0 * sind(theta_blade/2);
  Rl = (dl^2 - r0^2 + (c/2 + r0 * cosd(theta_blade/2))^2) / (2 * (dl + r0));
  ycl = y0 - tb/2 - Rl;

  xl = x;
  yl = ycl + sqrt(Rl^2 - x.^2);

  %smooth edges
  y_edge0_dx = r0 * sind(theta_blade/2);
  x_edge0_dx = + (c/2 - r0 * cosd(theta_blade/2));

  theta_r0 = 0:360;%-100:65;
  x_edge_dx = r0 * cosd(theta_r0) + x_edge0_dx;
  y_edge_dx = r0 * sind(theta_r0) + y_edge0_dx;

  y_edge0_sx = r0 * sind(theta_blade/2);
  x_edge0_sx = - (c/2 - r0 * cosd(theta_blade/2));

  theta_r0 = 0:360;%115:260;
  x_edge_sx = r0 * cosd(theta_r0) + x_edge0_sx;
  y_edge_sx = r0 * sind(theta_r0) + y_edge0_sx;

  
  %% Output vector
  x_profile = [flip(xu), xl]';
  y_profile = [flip(yu), yl]';
  
  % Flip if theta < 0
  if theta < 0
    
    y_profile = y_profile - 2 * y_profile;
    yl = yl - 2 * yl;
    yu = yu - 2 * yu;
    y = y - 2 * y;
    
    
    
  end
  
  % Print Blades
  figure
  plot(x,y,'k--')
  hold on
  plot(xu,yu,'k','linewidth',2)
  plot(xl,yl,'k','linewidth',2)
  %plot(x_edge_dx,y_edge_dx)
  %plot(x_edge_sx,y_edge_sx)
  title(['DCA Profile - c = ', num2str(c),', t = ',num2str(tb),'c'])
  axis equal



else
  return

  
end  
 


end
