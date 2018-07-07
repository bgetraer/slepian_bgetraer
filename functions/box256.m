function [ lon2N, lat2N ] = box2N(lond,latd,n)
%BOX2N Resample a square scattered, possible warped grid on a 2^nx2^n square
%mesh.

refvector = linspace(1,size(lond,1),2^n);

[refx, refy] = ndgrid(1:size(lond,1),1:size(lond,2));
[refx2N,refy2N] = ndgrid(refvector,refvector);

lon2N = griddata(refx,refy,lond,refx2N,refy2N);
lat2N = griddata(refx,refy,latd,refx2N,refy2N);
end

