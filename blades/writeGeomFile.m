function writeGeomFile(x,y,z,name)

% WRITEGEOMFILE writes a .dat file containing (x,y,z) coordinates
%
% This is meant to create a ready-to-use .txt for CAD import
%
% Example:
%
%   writeGeomFile(x,y,z,naca0012)

name = strcat(name,'.txt');

strPath = [pwd '/blades/points'];
strFull = fullfile(strPath,name);

n = length(x);
npoints = num2str(n);

if isempty(z)
    z = zeros(length(x),1);
elseif isempty(y)
    y = zeros(length(x),1);
elseif isempty(x)
    x = zeros(length(z),1);
end


fid = fopen(strFull,'wt');
%fprintf(fid, npoints);

for k = 1:n
    
    fprintf(fid, '%f   %f   %f\n',x(k), y(k), z(k));
    
end

fclose(fid);

end