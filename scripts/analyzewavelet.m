%**************************************************************************
%   Script developed for "REGIONAL FORCING OF GREENLAND ICE LOSS 2002-2017"
%   Spring 2018 Junior Paper, Princeton Department of Geosciences
%
%   SCRIPT 4
%   Takes a chosen wavelet decomposition and eliminates coefficients
%   based on contribution to a test image (here, the TOTAL change in 
%   Greenland mass over the GRACE time series). In other words, if the
%   wavelet doesn't represent anything that is noticeably changing,
%   eliminate it. This script thresholds wavelet coefficients by magnitude,
%   requiring the reconstruction to maintain 90% pixel-to-pixel invariance 
%   (as defined by IMINVAR) from the original image.
%   Saves an index for wavelet coefficients which pass the thresholding.
%   
%   PREVIOUS: CHOOSEWAVELET.m
%   NEXT: WAVEINPOLY.m
%
%   Benjamin Getraer bgetraer@princeton.edu
%   Modified: 2/21/2019
%   See also: PRCTILETHOLD, IMINVAR, WAVEDEC2, WAVEREC2, WTHCOEF2,
%           HISTWTHCOEF
%**************************************************************************

% locate slepian_bgetraer function and datafile directories, and set workspace
homedir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/';
functiondir = fullfile(homedir,'functions');
datadir = fullfile(homedir,'datafiles');
addpath(functiondir,datadir);   clear('homedir','functiondir');
setworkspace();


% load datafiles created in BOXGREENLAND.m and IMAGERYSEQ.m
load ptsGL 
load im_tools 
load im_seqSH
%% THRESHOLD WAVELETS ON THE FULL TIMESERIES DIFFERENCE IMAGE BY MAGNITUDE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Caluculate coefficient percentile threshold, create an index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Our reference image is the difference between the first and last date in 
%   the timeseries
refIm =  D(:,:,end)-D(:,:,1);
wavename = 'haar';
invarT = 0.9;                  % invariance target

% Index for the wavelet coefficients of our reference image which pass the 
% threshold determined by the invariance target provided.
[~,~, threshpassindex] = prctileThold(refIm, invarT, wavename);
% save for use in other files
save(fullfile(datadir,'threshpassindex'),'threshpassindex')
%% PLOTTED EXAMPLES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Visualize what we have done
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target = [0.90 0.99 1]; % the invariance targets
% originalimage = imread('lena_std.tif'); % the Lena test image
% originalimage = double(originalimage(:,:,3));
originalimage = refIm; % our Greenland difference map

figure(1)
clf

for i = 1:3
    subplot(2,3,i)
    [recon_image, prctGone] =...
        histWTHCOEF(originalimage,target(i),wavename);
    axis tight
    set(gca,'xlim',[-15,12])
    plottitle = sprintf('~%i%% invariance',round(target(i)*100));
    title(plottitle)
    
    subplot(2,3,i+3)
    imshow(recon_image,[])
    colormap(gca,bluewhitered([],1))
    plottitle = sprintf('%0.2f%% coefficients eliminated',round(prctGone*100,2));
    title(plottitle)
    hold on
    plot(gx,gy,'k','linewidth',1)
end