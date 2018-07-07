function [ lon2N, lat2N ] = box2N(lond,latd,n)
%BOX2N Resample a square scattered, possible warped grid on a 2^nx2^n square
%mesh.
%
% INPUT
%   lond    generally, the first dimension of the original grid (across)
%   latd    generally, the second dimention of the original grid (up)
%   n       the side lendgth of the grid is 2^n -- default is 256x256
%
% OUTPUT
%   lon2N   the resampled first dimension of the grid
%   lat2N   the resampled second dimension of the grid
%
% See Also: CUBE2SPHERE, BOX2N
%
% last modified by: bgetraer@princeton.edu, 6/25/2018


defval('n',8);  % default is a 256x256 box

refvector = linspace(1,size(lond,1),2^n);

[refx, refy] = ndgrid(1:size(lond,1),1:size(lond,2));
[refx2N,refy2N] = ndgrid(refvector,refvector);

lon2N = griddata(refx,refy,lond,refx2N,refy2N);
lat2N = griddata(refx,refy,latd,refx2N,refy2N);
end

