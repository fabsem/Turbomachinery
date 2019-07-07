function [x,y] = profilePoints(profileName,c,camberAngle,n)

%PROFILEPOINTS Evaluate Profile Points for Compressor blades
%
% Output vectors order:
%
%   x = [TE to LE upper side, LE to TE lower side];
%   y = [TE to LE upper side, LE to TE lower side];
%
% Input Data Required:
%
% profileName: Select the profile between 'naca65' and 'MCA'
% c: profile chord
% camberAngle: Profile Camber angle
% n: Number of points
%
% Example_
%
%   [x,y] = profilePoints('naca65',1,20,16)

  %Origin for radius
  x0 = 0;
  y0 = 0;

  %Chord
  %c = 1;

  if strcmp(profileName,'naca65')

      %NACA 65 Points in %chord
      yt = [0 1.124 1.571 2.222 2.709 3.111 3.746 4.218 4.824 5.057 4.87 4.151...
            3.038 1.878 0.749 0.354 0.15]/100;

  elseif strcmp(profileName,'MCA')

            %MCA Points
    yt = [0.1 0.217 0.334 0.534 0.778 0.985 1.369 1.711 2.275 2.678...
    2.9395 3 2.819 2.274 1.364 0.771 0.1]/100;

  else

    disp('No profile entered or non available series. Switching to NACA 65 as default')
    yt = [0 1.124 1.571 2.222 2.709 3.111 3.746 4.218 4.824 5.057 4.87 4.151...
    3.038 1.878 0.749 0.354 0.15]/100;

  end

  %Camber line slope
  phic = deg2rad(camberAngle);

  % Radius of mean camber line
  Rc = c / (2 * sin(phic));

  % SISTEMA: FAI CHE PUOI CAMBIARE N
  if n ~= 17
    % Number of nodes
    n = 17;

  end

  % Order of the node
  j = 1:n;

  % Mean Camber line coordinates
  for i = 1:length(j)
    
    xc(1,i) = x0 + j(i) * c / (n - 1);
    yc(1,i) = y0 + sqrt(Rc^2 - ( xc(i) - x0 )^2);

  end

  % Profile Coordinates
  % Upper part
  xu = xc - yt * sin(phic);
  yu = yc + yt * cos(phic);
  % Lower part
  xl = xc + yt * sin(phic);
  yl = yc - yt * cos(phic);

  % Rotation
  rot = pi*1/8;
  R = [cos(rot) -sin(rot); sin(rot) cos(rot)];
  c_coord = R * [xc; yc];
  l_coord = R * [xl; yl];
  u_coord = R * [xu; yu];

  % Rotated coordinates
  xc = c_coord(1,:);
  yc = c_coord(2,:);

  xl = l_coord(1,:);
  yl = l_coord(2,:);

  xu = u_coord(1,:);
  yu = u_coord(2,:);

  plot(xc,yc)
  hold on
  plot(xu,yu)

  plot(xl,yl)
  hold off
  axis equal

  %output
  x = [flip(xu)'; (xl)'];
  y = [flip(yu)'; (yl)'];

end
