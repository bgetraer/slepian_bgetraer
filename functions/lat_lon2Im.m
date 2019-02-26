function [ Fx,Fy,xp,yp ] = lat_lon2Im(lond,latd)
%LAT_LON2IM Create interpolant functions that transform global lat/lon
%   coordinates into their image basis.
%
% INPUT
%   lond     the longitude matrix
%   latd     the latitude matrix
%
% OUTPUT
%   Fx      the X dimension interpolant function (lat/lon to X)
%   Fy      the Y dimension interpolant function (lat/lon to Y)
%
% See Also: SCATTEREDINTERPOLANT, NDGRID
%
% last modified by: bgetraer@princeton.edu, 6/25/2018

% create a temporary image to extract properties
tempfig = figure(4);
clf
im1 = imagesc(lond); % the image representation of the matrix
% the number of x and y points
xpts = im1.XData(2);
ypts = im1.YData(2);
% the limits of the image axis
xylim = axis;
close(tempfig)

% the grid of x and y coordinates of the image basis
xp = linspace(xylim(1),xylim(2),xpts);
yp = linspace(xylim(3),xylim(4),ypts);
[xp,yp] = ndgrid(xp,yp);
xp=xp';
yp=yp';

% Interpolant functions for (lon,lat) to x and y
Fx = scatteredInterpolant(lond(:),latd(:),xp(:));
Fy = scatteredInterpolant(lond(:),latd(:),yp(:));
end

