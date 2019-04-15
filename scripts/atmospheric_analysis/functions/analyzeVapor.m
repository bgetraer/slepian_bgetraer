%**************************************************************************
% Script for analyzing Water Vapor over the Greenland Ice Sheet
%
% Last modified by: bgetraer@princeton.edu 3/30/2019
%**************************************************************************
% locate slepian_bgetraer function and datafile directories, and set workspace
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
datadir = fullfile(dir,'datafiles/');
merradir = fullfile(datadir,'MERRA2');
matDir = fullfile(merradir,'MerraMat');
addpath(dir,datadir);
addpath(fullfile(dir,'scripts/atmospheric_analysis/functions'))
addpath(fullfile(dir,'functions'))
setworkspace();

QV = load(fullfile(matDir,'Vapor2003-2010'),'QVdata','t','spacelim');
QV.QVdata = squeeze(QV.QVdata);
%%
startdate = datenum(2003,01,01);
enddate = datenum(2010,12,31);
alldates = startdate:enddate;

QVdata500 = squeeze(QV.QVdata(:,:,3,:));

figure(2)
clf
[QV,meanQV,QVA,fp,fl,ml,x,y] = plotoverGRIS(alldates,QVdata500 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,1)
hold on
ylabel('kg/kg')
title('mean 500 hPa QV over GRIS')

subplot(3,3,4:5)
hold on
ylabel('kg/kg')
title('mean 500 hPa QV anomaly over GRIS')

%% NORMALIZE PIXEL BY PIXEL
[ QVslopes ] = linearMap( QVdata500, alldates );

%% PLOT
figure(5)
clf
imj = QVslopes;
imPlot(imj.*365)
colorbar
axis square off
