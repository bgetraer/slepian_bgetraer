%**************************************************************************
%   JP02    Evaluate slepian solutions for Greenland on a grid using the
%   cubed sphere.
%
%
%   See also: PLM2XYZ CUBE2SPHERE
%**************************************************************************

addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');
%% In an order 60 Spherical harmonic system,
%   the highest spatial period is 6 degrees lon, 3 degrees lat.
% The nyquist sampling is every 3 degrees lon, every 1.5 degrees lat.
% Greenland is ~60-85 deg lat, 20-70 deg lon. or 50x25 degrees.
%   50/3 = 16.667 25/1.5 = 16.667
%
% x = 0:1/(60*2/180):180;
% y = cos(60*2/180*pi*x);
% x2 = 0:0.001:180;
% y2 = cos(60*2*pi/180*x2);
% clf
% hold on
% plot(x2,y2)
% plot(x,y,'x')
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
order = 11;  % # points on cubed sphere = 6*2^order
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
%% get the Greenland data in lmcosi matrix
L=60;

load('Greenland60data');
N=19.6735;
signal = G(:,1:20)*(slepcoffs(1,1:20)-mean(slepcoffs(:,1:20)))';

% Create blank LMCOSI matrix
[~,~,~,blank,~,~,~,~,~,ronm]=addmon(L);
lmcosi_mat = zeros([size(blank) size(signal,3)]) + blank;


% Create the coefficient blanks
cosi=blank(:,3:4);
% grab the coefficients of an alpha eigentaper and
% re-index them in lmcosi format
cosi(ronm)=signal(:,1,1);
% Add them to the full matrix
lmcosi_mat(:,3:4,1)=cosi;

% Create 2d matrix of summed coefficients for all grabbed alphas.
lmcosi_sum = [blank(:,1:2) sum(lmcosi_mat(:,3:4,:),3)];

%% CALCULATE AT SOME POINTS
% this is the time consuming bit... plm2xyz takes a while
[data,lon,lat,Plm]=plm2xyz(lmcosi_sum,latd(:),lond(:)); % vector of solutions
D = reshape(data,shp);  % put the vector back into matrix form
Dim = imrotate(D,180);
%%  MORE PLOTTING
% load(filename)
% 
% figure(3)
% clf
% plotcont([],[],3)
% hold on
% [s1,s2,s3]= sphere;
% load topo
% surf(s1,s2,s3,'edgecolor','none','FaceColor',[.9 .9 .9])
% surface(xprime,yprime,zprime,'edgecolor','none','FaceColor','texture','Cdata',D);
% axis()
% colormap(bluewhitered(1000))
% colorbar
% 
% figure(4)
% clf
% imagesc(Dim)
% axis image
% colormap(bluewhitered(1000))
% colorbar
%% SAVE DATA
% filename1 = sprintf('boxGL%d',order);
% save(fullfile(datadir,filename1),'D','Dim','xprime','yprime','zprime','latd','lond')
filename2 = sprintf('ptsGL%d',order);
save(fullfile(datadir,filename2),'xprime','yprime','zprime','latd','lond')
