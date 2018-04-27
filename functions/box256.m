function [ lon256, lat256 ] = box256(lond,latd)
%BOX256 Resample a square scattered, possible warped grid on a 256x256 square
%mesh.

refvector = linspace(1,size(lond,1),256);

[refx, refy] = ndgrid(1:size(lond,1),1:size(lond,2));
[refx256,refy256] = ndgrid(refvector,refvector);

lon256 = griddata(refx,refy,lond,refx256,refy256);
lat256 = griddata(refx,refy,latd,refx256,refy256);

end

