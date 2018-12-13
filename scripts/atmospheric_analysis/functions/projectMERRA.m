function [ Flon2x,Flat2y,gx,gy,bx,by,contx,conty ] = projectMERRA( Vdata, thespacelim )
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

% x and y limits in image basis
xIMlim = [0 size(Vdata,1)]+0.5;
yIMlim = [0 size(Vdata,2)]+0.5;

% Interpolant Functions
Flon2x = griddedInterpolant(thespacelim(1,:),xIMlim);
Flat2y = griddedInterpolant(thespacelim(2,:),flip(yIMlim));


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
end

