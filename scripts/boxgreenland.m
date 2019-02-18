%**************************************************************************
%   Script developed for "REGIONAL FORCING OF GREENLAND ICE LOSS 2002-2017"
%   Spring 2018 Junior Paper, Princeton Department of Geosciences
%
%   Includes script for FIGURE 5
%  
%   SCRIPT 1
%   Create grid for Greenland on a square 2^n grid using the cubed sphere.
%   This is the first script in a series that all accomplish small
%   pieces of the puzzle. 
%   NEXT: IMAGERYSEQ.m
%
%   Benjamin Getraer bgetraer@princeton.edu
%   Modified: 6/25/2018
%   See also: PLM2XYZ, CUBE2SPHERE, SQUAREGRID, BOX2N, PLOTPLMGREENLAND,
%   LAT_LON2IM, JP02FIG5
%**************************************************************************

% locate slepian_bgetraer function and datafile directories, and set workspace
homedir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/';
functiondir = fullfile(homedir,'functions');
datadir = fullfile(homedir,'datafiles');
addpath(functiondir,datadir);   clear('homedir','functiondir');
setworkspace();

%% SET UP THE GRID AROUND GREENLAND

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Create the square grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the Euler angles to center around Greenland
alpha = 0;
beta = -0.31;
gamma = 0.79;

% Create a 256x256 grid around Greenland
%   returns the coordinates of the grid points in lat/lon and cartesian
[lond, latd, xprime, yprime, zprime] = squareGrid([],alpha,beta,gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global lat/lon to image basis interpolant functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ Fx,Fy ] = lat_lon2Im(lond,latd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project specific Lat/Lon onto the Image basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Greenland lat/lon coordinates
gxy = greenland(10);
buff = greenland(10,0.5);
% Greenland area transformation
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
jp02fig5(xprime,yprime,zprime,lond,latd,gx,gy,bx,by,Fx,Fy)

%% SAVE DATA TO IMPORT INTO OTHER SCRIPTS
% The global grid points around Greenland
filename1 = 'ptsGL.mat';
save(fullfile(datadir,filename1),'xprime','yprime','zprime','latd','lond')

% The image grid and tools around Greenland
filename2 = 'im_tools.mat';
save(fullfile(datadir,filename2),'Fx','Fy','gx','gy','bx','by','areaGL')

%**************************************************************************
%   End of BOXGREENLAND
%**************************************************************************
