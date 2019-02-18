function [handle] = flowArrow(u,v,ds,scale,linespec)
%FLOWARROW 
%
% INPUT
%   u   east component (straight from the MERRA)
%   v   north component (straight from the MERRA)
%   ds  downsample percent

xd = linspace(1,size(u,1),size(imresize(u,ds),1));
yd = linspace(1,size(v,2),size(imresize(v,ds),2));

u = imrotate(u,90);
v = imrotate(v,90);
handle = quiver(xd,yd,imresize(u,ds),-imresize(v,ds),scale,linespec);
end

