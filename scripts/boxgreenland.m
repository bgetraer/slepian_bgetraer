%**************************************************************************
%   JP02    Evaluate slepian solutions for Greenland on a grid using the
%   cubed sphere.
%
%   See also: PLM2XYZ CUBE2SPHERE
%**************************************************************************

addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

%% Set up the cubed sphere around Greenland
% greenland polygon, should we want it
[gxy] = greenland();
geom = polygeom(gxy(:,1),gxy(:,2));
% the centroid, should we want it
c = geom(2:3);

% convert from lat-lon to XYZ
gxy = gxy.*(2*pi/360);
c = c.*(2*pi/360);
[gx,gy,gz] = sph2cart(gxy(:,1),gxy(:,2),ones(size(gxy(:,1))));
[cx,cy,cz] = sph2cart(c(1),c(2),1);

% Euler angles to center around Greenland
alpha = 0;
beta = -0.31;
gamma = 0.79;
order = 10;  % # points on cubed sphere = 6*2^order
% Get the cubed sphere points
[x,y,z] = cube2sphere(order,alpha,beta,gamma,0,0);
% index the points that are around greenland
l = length(x);
indx = 11*l/32:21*l/32;
% indy = 13*l/32:19*l/32;
indy = 11*l/32:21*l/32;
% define the points that are around greenland
xprime = x(indx,indy,5);
yprime = y(indx,indy,5);
zprime = z(indx,indy,5);

% get the points in degrees lat lon
[lon,lat,r]=cart2sph(xprime,yprime,zprime); % the points in radians
shp = size(lon);            % the shape of the grid matrix
latd = lat*360/(2*pi);      % vector of lat in degrees
lond = 360+lon*360/(2*pi);  % vector of lon in degrees
% flip it around to be "right-side" up
lond = flip(flip(lond),2);
latd = flip(flip(latd),2);

% RESAMPLE TO GET A 256x256 square GRID ON GREENLAND!!
[lond, latd] = box256(lond,latd);
r = ones(256);
[xprime,yprime,zprime] = sph2cart(deg2rad(lond), deg2rad(latd),r);

figure
hold on
% plot(lond,latd,'k.')
% plot(lond, latd,'r.')

% Plot the points
clear p
figure(2)
clf
p=plot3(xprime(:),yprime(:),zprime(:),'o',...
    'MarkerF','r','MarkerE','r');
hold on; 
% set(p,'MarkerS',1)
plot3(xprime(1,end),yprime(1,end),zprime(1,end),'o',...
    'MarkerF','r','MarkerE','k');
plot3(xprime(1,1),yprime(1,1),zprime(1,1),'o',...
    'MarkerF','r','MarkerE','k');
plot3(xprime(end,1),yprime(end,1),zprime(end,1),'o',...
    'MarkerF','r','MarkerE','k');
plot3(xprime(end,end),yprime(end,end),zprime(end,end),'o',...
    'MarkerF','r','MarkerE','k');
plot3(gx,gy,gz,'linewidth',5);
plot3(cx,cy,cz,'x');
axis equal; hold off; axis off;
view(37,90)

%% SAVE DATA
% filename1 = sprintf('boxGL%d',order);
% save(fullfile(datadir,filename1),'D','Dim','xprime','yprime','zprime','latd','lond')
filename2 = sprintf('ptsGL%d',order);
save(fullfile(datadir,filename2),'xprime','yprime','zprime','latd','lond')


load(fullfile(datadir,sprintf('ptsGL%d.mat',order)))

% plot greenland and the outline of the box
figure(3)
clf 
hold on
greenland
% the outline
plot(lond(:,1),latd(:,1),'k-');
plot(lond(:,end),latd(:,end),'k-');
plot(lond(1,:),latd(1,:),'k-');
plot(lond(end,:),latd(end,:),'k-');
% the endpoints
plot(lond(1,1),latd(1,1),'o',...
    'MarkerF','r','MarkerE','k');
plot(lond(end,end),latd(end,end),'o',...
    'MarkerF','r','MarkerE','k');
plot(lond(1,end),latd(1,end),'o',...
    'MarkerF','r','MarkerE','k');
plot(lond(end,1),latd(end,1),'o',...
    'MarkerF','r','MarkerE','k');

%% Image of Greenland mass and image basis
figure(4)
clf
im1 = imagesc(lond); % the image of Greenland
hold on
axis image
colormap(bluewhitered(1000,0));
colorbar

% the number of x and y points
xpts = im1.XData(2);
ypts = im1.YData(2);

% the limits of the image axis
xylim = axis;

% the grid of x and y coordinates of the image basis
xp = linspace(xylim(1),xylim(2),xpts);
yp = linspace(xylim(3),xylim(4),ypts);
[xp,yp] = ndgrid(xp,yp);
xp=xp';
yp=yp';

%% Interpolant functions for (lon,lat) to x and y
Fx = scatteredInterpolant(lond(:),latd(:),xp(:));
Fy = scatteredInterpolant(lond(:),latd(:),yp(:));
%% Project Lat/Lon onto the Image
% Greenland lat lon
gxy = greenland(10);
buff=greenland(10,0.5);
% Greenland xy in the image basis
gx = Fx(gxy(:,1),gxy(:,2));
gy = Fy(gxy(:,1),gxy(:,2));
bx = Fx(buff(:,1),buff(:,2));
by = Fy(buff(:,1),buff(:,2));
% endpoints
x1=Fx(lond(1,1),latd(1,1));
y1=Fy(lond(1,1),latd(1,1));
x2=Fx(lond(1,end),latd(1,end));
y2=Fy(lond(1,end),latd(1,end));
x3=Fx(lond(end,1),latd(end,1));
y3=Fy(lond(end,1),latd(end,1));
x4=Fx(lond(end,end),latd(end,end));
y4=Fy(lond(end,end),latd(end,end));


figure(4)
clf
im1 = imagesc(lond);
axis image

colormap(bluewhitered(1000,0));
colorbar
hold on
% endpoints
plot([x1,x2,x3,x4],[y1,y2,y3,y4],'o',...
    'MarkerF','r','MarkerE','k');
% Greenland
plot(gx,gy,'k-')

%% save interpolant functions
save(fullfile(datadir,sprintf('im_tools%d.mat',order)),'Fx','Fy','gx','gy','bx','by','xp','yp')


%% ************************************************************************
%   ABANDONED PLOTTING ROUTINES
% % *************************************************************************
% %% In an order 60 Spherical harmonic system,
% %   the highest spatial period is 6 degrees lon, 3 degrees lat.
% % The nyquist sampling is every 3 degrees lon, every 1.5 degrees lat.
% % Greenland is ~60-85 deg lat, 20-70 deg lon. or 50x25 degrees.
% %   50/3 = 16.667 25/1.5 = 16.667
%
% x = 0:1/(60*2/180):180;
% y = cos(60*2/180*pi*x);
% x2 = 0:0.001:180;
% y2 = cos(60*2*pi/180*x2);
% clf
% hold on
% plot(x2,y2)
% plot(x,y,'x')
%
% %% get the Greenland data in lmcosi matrix
% L=60;
% 
% load('Greenland60data');
% N=19.6735;
% signal = G(:,1:20)*(slepcoffs(1,1:20)-mean(slepcoffs(:,1:20)))';
% 
% % Create blank LMCOSI matrix
% [~,~,~,blank,~,~,~,~,~,ronm]=addmon(L);
% lmcosi_mat = zeros([size(blank) size(signal,3)]) + blank;
% 
% 
% % Create the coefficient blanks
% cosi=blank(:,3:4);
% % grab the coefficients of an alpha eigentaper and
% % re-index them in lmcosi format
% cosi(ronm)=signal(:,1,1);
% % Add them to the full matrix
% lmcosi_mat(:,3:4,1)=cosi;
% 
% % Create 2d matrix of summed coefficients for all grabbed alphas.
% lmcosi_sum = [blank(:,1:2) sum(lmcosi_mat(:,3:4,:),3)];
% 
% %% CALCULATE AT SOME POINTS
% % this is the time consuming bit... plm2xyz takes a while
% [data]=plm2xyz(lmcosi_sum,latd(:),lond(:)); % vector of solutions
% D = reshape(data,shp);  % put the vector back into matrix form
% % Dim = imrotate(D,180);
% %%  MORE PLOTTING
% % load(filename)
% % 
% % figure(3)
% % clf
% % plotcont([],[],3)
% % hold on
% % [s1,s2,s3]= sphere;
% % load topo
% % surf(s1,s2,s3,'edgecolor','none','FaceColor',[.9 .9 .9])
% % surface(xprime,yprime,zprime,'edgecolor','none','FaceColor','texture','Cdata',D);
% % axis()
% % colormap(bluewhitered(1000))
% % colorbar
% % 
% % figure(4)
% % clf
% % imagesc(Dim)
% % axis image
% % colormap(bluewhitered(1000))
% % colorbar