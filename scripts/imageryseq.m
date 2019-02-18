%**************************************************************************
%   Script developed for "REGIONAL FORCING OF GREENLAND ICE LOSS 2002-2017"
%   Spring 2018 Junior Paper, Princeton Department of Geosciences
%
%   SCRIPT 2
%   Evaluate the GRACE spherical harmonic functions on the image grid for 
%   each year in the time-series. ** WARNING: PLM2XYZ takes a lot of time
%   and memory while running. Do not run the entire time-series unless you
%   need it. For an example of how the script and images work, choose only 
%   a few of the dates. ** A plotting routine of an existing time-series
%   datafile can be run independently (see the end of the script).
%   This is the second script in a series that all accomplish small
%   pieces of the puzzle. 
%   PREVIOUS: BOXGREENLAND.m
%   NEXT: CHOOSEWAVELET.m
%
%   Benjamin Getraer bgetraer@princeton.edu
%   Modified: 7/5/2018
%   See also: BOXGREENLAND.m, PLM2GRID, GRACE2PLMT, PLM2XYZ
%**************************************************************************

% locate slepian_bgetraer function and datafile directories, and set workspace
homedir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/';
functiondir = fullfile(homedir,'functions');
datadir = fullfile(homedir,'datafiles');
addpath(functiondir,datadir);   clear('homedir','functiondir');
setworkspace();

% load datafiles created in BOXGREENLAND.m
load('ptsGL')
load('im_tools')

%% IMPORT THE GRACE DATA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Import the GRACE surface density coefficients, 2003--2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the CSR RL05 L=60 GRACE surface density coefficients
[sdcoffs,~,thedates]=grace2plmt('CSR','RL05','SD',0);

% crop off first 7 months of untrustworthy data (see Harig & Simons, 2016): 
validrange = monthnum(1,2003,thedates):length(thedates); 
% define the surface density spherical harmonic coefficients and dates
sdcoffs = sdcoffs(validrange,:,:);
thedates = thedates(validrange);

%% Evaluate the signals on the grid
%   this is time consuming... plm2xyz takes a while
%   data is saved as datafile im_seqSH.mat
D = plm2grid(sdcoffs,thedates,latd,lond);
%% Plot the sequence as a movie
% Load an existing file, or the one you just made
load(fullfile(datadir,'im_seqSH'))
% The differences between an image and the first image (Jan 2003)
Dprime = D-D(:,:,1);
% Plot as a movie
figure()
for i=1:size(D,3)
    clf
    imagesc(Dprime(:,:,i));
    caxis([-3.8344    1.5619]*1E3)
    colormap(bluewhitered(1000,1));
    colorbar
    axis image
    hold on;axis off;
    % Greenland
    plot(gx,gy,'k-')
    title(sprintf('%s',datestr(thedates(i),'mmm YYYY')))
    pause(0.1)
end

%**************************************************************************
%   End of IMAGERYSEQ
%**************************************************************************