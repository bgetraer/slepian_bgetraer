addpath /Users/benjamingetraer/Downloads/
DEM  = imread('GIMP_dem.tif');
geotiffread
%%
DEMsmall = imresize(DEM,1/3);
%%
% x and y limits in image basis
xIMlim = [0 size(DEMsmall,2)]+0.5;
yIMlim = [0 size(DEMsmall,1)]+0.5;

thespacelim = [-75, -14;60,83];
% Interpolant Functions
Fx = griddedInterpolant(thespacelim(1,:),xIMlim);
Fy = griddedInterpolant(thespacelim(2,:),flip(yIMlim));

gxy = greenland(10);

gx = Fx(gxy(:,1)-360);
gy = Fy(gxy(:,2));
%%
figure(1)
clf
imshow(DEMsmall,[])
hold on
plot(gx,gy,'r-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
[DEM,R]  = geotiffread('GIMP_dem.tif');

%%
DEMsmall = imresize(DEM,1/3);

%%
% x and y limits in image basis
xIMlim = [0 size(DEMsmall,2)]+0.5;
yIMlim = [0 size(DEMsmall,1)]+0.5;

thespacelim = [-75, -14;60,83];
% Interpolant Functions
Fx = griddedInterpolant(thespacelim(1,:),xIMlim);
Fy = griddedInterpolant(thespacelim(2,:),flip(yIMlim));

gxy = greenland(10);

[gx,gy] = worldToIntrinsic(R,gxy(:,1),gxy(:,2));
%%
figure(1)
clf
imshow(DEM,[])
%%
hold on
% plot(gx,gy,'r-')
origin = [-55.7983, 59.1996];
plot(worldToIntrinsic(R,0,0),'*r')