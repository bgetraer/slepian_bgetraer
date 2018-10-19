%**************************************************************************
%   Script altered for SENIOR THESIS, Princeton Department of Geosciences
%
%  
%   SCRIPT 1.b
%   Create grid for Greenland on a square 2^n grid using the cubed sphere.
%   This is the first script in a series that all accomplish small
%   pieces of the puzzle. THIS SCRIPT IS BOXGREENLAND WITH THE BASIS
%   ROTATED ~90 DEGREES AND SIZED 128x128
%   NEXT: IMAGERYSEQ_B.m
%
%   Benjamin Getraer bgetraer@princeton.edu
%   Modified: 6/25/2018
%   See also: BOXGREENLAND, PLM2XYZ, CUBE2SPHERE, SQUAREGRID, BOX2N, 
%   PLOTPLMGREENLAND, LAT_LON2IM, JP02FIG5
%**************************************************************************

addpath('/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/IndependentWork/SH_Workspace');

%% SET UP THE GRID AROUND GREENLAND

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Create the square grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the Euler angles to center around Greenland
alpha = 0.8;
beta = -0.35;
gamma = 0.75;


% Create a 128x128 grid around Greenland
%   returns the coordinates of the grid points in lat/lon and cartesian
n=7; % image with 2^n points
[lond, latd, xprime, yprime, zprime] = squareGrid([],alpha,beta,gamma,[],[],n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global lat/lon to image basis interpolant functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ Fx,Fy ] = lat_lon2Im(lond,latd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project specific Lat/Lon onto the Image basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Greenland lat/lon coordinates
gxy = greenland(10);
buff=greenland(10,0.5);
% Greenlan area transformation
% Since Int should have units of (fn * m^2), need to go from fractional
%   sphere area to real area.  If the fn is surface density, this output is
%   in kilograms.  Then change the units from kg to Gt in METRIC tons.
areaGL = areaint(buff(:,2),buff(:,1))*4*pi*fralmanac('a_EGM96','Earth')^2;

% Greenland xy in the image basis
gx = Fx(gxy(:,1),gxy(:,2));
gy = Fy(gxy(:,1),gxy(:,2));
bx = Fx(buff(:,1),buff(:,2));
by = Fy(buff(:,1),buff(:,2));

%% FIGURE: GLOBAL v IMAGE BASIS:
jp02fig5(xprime,yprime,zprime,lond,latd,gx,gy,bx,by,Fx,Fy,n)

%% SAVE DATA TO IMPORT INTO OTHER SCRIPTS
% The global grid points around Greenland
filename1 = 'ptsGL_b.mat';
save(fullfile(datadir,filename1),'xprime','yprime','zprime','latd','lond')

% The image grid and tools around Greenland
filename2 = 'im_tools_b.mat';
save(fullfile(datadir,filename2),'Fx','Fy','gx','gy','bx','by','areaGL')