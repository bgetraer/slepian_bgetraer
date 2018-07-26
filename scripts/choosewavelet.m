%**************************************************************************
%   Script developed for "REGIONAL FORCING OF GREENLAND ICE LOSS 2002-2017"
%   Spring 2018 Junior Paper, Princeton Department of Geosciences
%  
%   Includes script for FIGURE 6
%
%   SCRIPT 3
%   Tests various 2 dimensional wavelet decompositions on a test image in
%   order to compare relative tradeoffs in variance and bias after partial
%   reconstruction. The test image chosen in my example is the TOTAL change
%   in Greenland mass over the GRACE time series.
%   PREVIOUS: IMAGERYSEQ.m
%   NEXT: ANALYZE WAVELET.m
%
%   Benjamin Getraer bgetraer@princeton.edu
%   Modified: 7/25/2018
%   See also: JP02FIG6, WAVEFUN, WAVEDEC2, WTHCOEF2, WAVEREC2
%**************************************************************************

addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

load('ptsGL')
load('im_tools')
load('im_seqSH')
%% Parameters for the wavelet testing

% Define test image (here the difference between the first date and last
%   date of the GRACE time series. "D" is a datafile produced in PLM2GRID
%   (see IMAGERYSEQ.m).
testimage = D(:,:,end)-D(:,:,1);

% Choose wavelets to test by MATLAB name. Here I chose to do 3 at a time
%   for clarity of plotting. In order, the Haar, the Fejer-Korovkin 4, and
%   the biorthoganal 3.7 wavelets are chosen for this example.
wname = {'haar','fk4','bior3.7'};
% Decomposition levels to include: 10 is the max for these wavelets on the
%   256x256 grid.
level = 10;

% Thresholding percentiles: what below what percentile of wavelet 
%   coefficient values do we throw away the data? Here I choose specific
%   ones that will plot nicely.
ptl = sort([linspace(96,100,50),99.8,97.906,99.745]);

% Do the decompositions and generate the figure
jp02fig6(testimage, wname, level, ptl)
