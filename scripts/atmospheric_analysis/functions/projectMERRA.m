function [ Flon2x,Flat2y,gx,gy,bx,by,contx,conty,azores,reykjavik, X,Y,LON,LAT] ...
    = projectMERRA( Vdata, thespacelim )
%PROJECTNETCDF Creates an interpolant function for projecting LAT LON onto
%the basis of the variable matrix image. 
%   Outputs the interpolant functions, and the locations of Greenland, a
%   0.5 deg buffer around Greenland, and the continents
% INPUT
%   Vdata       MERRA2 data matrix, first two dimensions are spatial
%   thespacelim the [lon;lat] limits of the matrix
% OUTPUT   
%   Flon2x      Interpolant functions: 
%   Flat2y          accepts lat/lon, returns image coordinates
%   gx,gy       coordinates for greenland
%   bx,by       coordinates for 0.5 deg buffer around greenland
%   contx,conty coordinates for continents
%   azores      coordinates for azores
%   reykjavik   coordinates for reykjavik
%
% SEE ALSO:
%   GRIDDEDINTERPOLANT
%
% Last modified by bgetraer@princeton.edu 3/9/2019


% x and y limits in image basis
xIMlim = [0 size(Vdata,1)]+0.5;
yIMlim = [0 size(Vdata,2)]+0.5;

% Interpolant Functions
Flon2x = griddedInterpolant(thespacelim(1,:),xIMlim);
Flat2y = griddedInterpolant(thespacelim(2,:),flip(yIMlim));

% Full meshgrid indexes
lonpoints = linspace(thespacelim(1),thespacelim(3),size(Vdata,1));
latpoints = linspace(thespacelim(2),thespacelim(4),size(Vdata,2));
xpoints = Flon2x(lonpoints);
ypoints = Flat2y(latpoints);
[X,Y]       = meshgrid(xpoints,ypoints);
[LON,LAT]   = meshgrid(lonpoints,latpoints);
LON(LON<0)  = LON+360;


% add GREENLAND and PLOTCONT function directory
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
addpath(fullfile(dir,'functions'))
setworkspace('/Users/benjamingetraer/Documents/IndependentWork/SH_Workspace');

% Get lat/lon coordinates
gxy = greenland(10);
buff= greenland(10,0.5);
[~,~,XYZ] = plotcont();
close

% Greenland xy in the image basis
gx = Flon2x(gxy(:,1)-360);
gy = Flat2y(gxy(:,2));
bx = Flon2x(buff(:,1)-360);
by = Flat2y(buff(:,2));
contx = Flon2x(XYZ(:,1)-360);
conty = Flat2y(XYZ(:,2));
azores = [Flon2x( -27.853785), Flat2y(38.471189)];
reykjavik = [Flon2x(-21.933333), Flat2y(64.133333)];
end

