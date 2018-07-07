function [lond, latd, xprime, yprime, zprime] = squareGrid(order,alpha,beta,gamma,face,size,n)
%SQUAREGRID A squareGrid for use in evaluating spherical harmonics for
%wavelet analysis.
%   Create a square grid around a given center using the cubed sphere as 
%   guidance.
%   The square grid can be rotated to any orientation around the center and 
%   can be of any resolution in powers of 2 with a default resolution of 
%   256x256.
%
% INPUT
%   order   # points on cubed sphere = 6*2^order
%               tradeoff in run-time vs resolution quality of square grid
%               following the cubed sphere grid lines. default is 10.
%   alpha, beta, gamma      Euler angles for the center of the square grid.
%                           These define the location and tilt of the grid.
%   face    the aligned face of the cubed sphere (1 through 6) default 5.
%   size    the size of the square around the center (not the size of the
%           grid!) relative to the cubed sphere. In other words, the number
%           of cubed square points that you created that you want the new
%           square grid side to span.  
%   n      The side length of the new square grid is 2^n
%                   (default 256x256, ie n=8)
% OUTPUT
%   lond, latd              the square grid in lat lon
%   xprime,yprime,zprime    the square grid in cartesian coordinates
%
% See Also: CUBE2SPHERE, BOX2N
%
% last modified by: bgetraer@princeton.edu, 6/6/2018

defval('order',10);  % # points on cubed sphere = 6*2^order
defval('size',320);  
defval('face',5);
defval('n',8);

% Get the cubed sphere points
[x,y,z] = cube2sphere(order,alpha,beta,gamma,0,0);
% the length of each side of the cubed sphere
l = 2^(order);
% we only want the points in some square surrounding the center
ind = (l/2-size/2):(l/2+size/2);
% define the points that are around our center
xprime = x(ind,ind,face);
yprime = y(ind,ind,face);
zprime = z(ind,ind,face);

% get the points in degrees lat lon
[lon,lat]=cart2sph(xprime,yprime,zprime); % the points in radians
latd = lat*360/(2*pi);      % vector of lat in degrees
lond = 360+lon*360/(2*pi);  % vector of lon in degrees
% flip it around to be "right-side" up
lond = flip(flip(lond),2);
latd = flip(flip(latd),2);

% RESAMPLE ALONG THE CUBED SPHERE GRID TO GET A NEW 2^N SQUARE GRID!!
[lond, latd] = box2N(lond,latd,n);  % the resampling function
r = ones(2^n);
[xprime,yprime,zprime] = sph2cart(deg2rad(lond), deg2rad(latd),r);

end

