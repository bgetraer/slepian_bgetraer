%**************************************************************************
%   Script developed for "REGIONAL FORCING OF GREENLAND ICE LOSS 2002-2017"
%   Spring 2018 Junior Paper, Princeton Department of Geosciences
%
%   Includes script for FIGURE 7
%
%   SCRIPT 5
%   Takes a chosen wavelet decomposition and eliminates coefficients
%   based on areal overlap with a polygon (GREENLAND). Saves the 
%   coefficient logical indexing and an image mask for reconstructing
%   images using only wavelets important within the polygon.
%   
%   PREVIOUS: ANALYZEWAVELET.m
%   NEXT: .m
%
%   Benjamin Getraer bgetraer@princeton.edu
%   Modified: 3/12/2019
%   See also: WAVELEVLINDEX, AREAINPOLYGON, JP02FIG7

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
load threshpassindex
%% Spatial index for points in and around Greenland

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MASK OF GREENLAND IN IMAGE BASIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz = size(xp);
a = inpolygon(xp(:),yp(:),bx,by);   % array of points inside Greenland
A = 1*reshape(a,sz);                % matrix of points inside Greenland

% plot a visualization
figure(1)
clf
imshow(A)
hold on
title("mask of image points inside of Greenland")
colormap(bone)
plot(bx,by,'r','linewidth',2)
plot(gx,gy,'k','linewidth',1)

save(fullfile(datadir,'GLimagemask'),'A')
%% Visualization of level indexing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   IMAGE RECONSTRUCTION BY DECOMPOSITION LEVEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
wavelevelINDEX('test');

%% Threshold by spatial support within Greenland

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   GENERATE UNIFORM WAVELET DECOMPOSITION DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavename = 'haar';
level = wmaxlev(sz,wavename);
[C,S]=wavedec2(A,level,wavename);
filename = 'passIndexbyArea.mat';

C0 = zeros(1,length(C));

if ~exist(fullfile(datadir,filename),'file')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   FIND THE WAVELETS WHICH PASS THE THRESHOLD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ I, L, coefinlevel] = wavelevelINDEX( C,S );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ polygonpassindex ] = ...
        areaInPolygon( C, S, wavename, level, xp, yp, bx, by, threshpassindex );
    
    save(fullfile(datadir,filename),'polygonpassindex')
else
    load(fullfile(datadir,filename))
end

%% SHOW SUCCESS OF MASK RECONSTRUCTION USING MAGNITUDE AND AREA THRESHOLDS

filename = 'passIndexbyAreaONLY';
load(fullfile(datadir,filename))

[C,S]=wavedec2(A,level,wavename);
PASS = waverec2(C.*polygonpassindexONLY,S,wavename);

figure(6)
clf
ax1 = subplot(1,2,1);
imagesc(PASS)
colormap(ax1,bone);
colorbar('westoutside')
axis image
hold on;
title('reconstruction of binary mask')
axis off

subplot(1,2,2)
imagesc(PASS-A)
colormap(jet);
colorbar('eastoutside')
axis image
hold on;
title('reconstruction error')
axis off

%% FIGURES OF TOTAL GRACE MELT, SH and WAVELET BASIS
GRACEdiff = D(:,:,end) - D(:,:,1);
[C,S]=wavedec2(GRACEdiff,level,wavename);
Cprime = C.*polygonpassindex;
GRACEdiffprime = waverec2(Cprime,S,wavename).*A;
GRACEdiffprime(~GRACEdiffprime)=nan;
figure(1)
clf
imagesc(GRACEdiff.*A,'AlphaData',~isnan(GRACEdiffprime));
hold on
axis image
plot(gx,gy,'k-')
colormap(bluewhitered([],1))
axis off

figure(2)
clf
imagesc(GRACEdiffprime,'AlphaData',~isnan(GRACEdiffprime));
hold on
axis image
plot(gx,gy,'k-')
colormap(bluewhitered([],1))
axis off

%% Wavelet reconstruction on the reference image using our thresholding

jp02fig7();

%% CHOOSE HAAR WAVELET, CHOOSE percentile threshold to minimize bias
wavename = 'haar';
target = 0.9;
figure(3)
clf

hold on

% Histogram plot:
[recon_image, prctGone] =...
        histWTHCOEF(GRACEdiff,target,wavename);
    axis tight
    set(gca,'xlim',[-15,12])
    plottitle = sprintf('~%i%% invariance',round(target*100));
    title(plottitle)