%**************************************************************************
% Script for comparing surface temperature and mass anomaly over the
% Greenland Icesheet
% Merra2 Data from:
% https://disc.gsfc.nasa.gov/datasets/M2I6NPANA_V5.12.4/summary?keywords=%22MERRA-2%22
%
%   SCRIPT 2
%
%   
%   PREV: ANALYZETEMP.m
%   NEXT:
%   
% SEE ALSO:
%   
%
% Last modified by: bgetraer@princeton.edu 3/13/2019
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

%% Load Temp data to a structure
meltData = load(fullfile(matDir,'normalizedMeltMap'));
% MERRA projection data
pM = load(fullfile(matDir,'projectMERRAGL'));
% GRACE projection data
pG = load('im_tools');
ptsGL = load('ptsGL');
% GRACE mass data
massData = load('im_seqSH');


% Get the MERRA LAT LON grid 
MerraLonX = pG.Fx(pM.LON,pM.LAT);
MerraLatY = pG.Fy(pM.LON,pM.LAT);
indexIN = inpolygon(MerraLonX,MerraLatY,pG.bx,pG.by);
% MerraLonX = MerraLonX(indexIN);
% MerraLatY = MerraLatY(indexIN);


GRACEdiff = massData.D(:,:,end) - massData.D(:,:,1);
mass = interp2(ptsGL.xp,ptsGL.yp,GRACEdiff,MerraLonX,MerraLatY)';
%%
figure(1)
clf
imagesc(GRACEdiff)
hold on
% plot(ptsGL.xp,ptsGL.yp,'kx')
plot(MerraLonX,MerraLatY,'kx')
plot(pG.gx,pG.gy,'rx')



%%
totalmelt = sum(meltData.allmeltMap,3);

figure(10)
clf
scatter(totalmelt(indexIN'),-mass(indexIN'),5,MerraLonX(indexIN'),'x')
colorbar

% figure(11)
% clf
% imagesc(MerraLonX(indexIN'))
% colorbar


%%

load im_seqSH
load passIndexbyArea
load GLimagemask


MerraLonX = Fx(LON,LAT);
MerraLatY = Fy(LON,LAT);

indexIN = inpolygon(MerraLonX,MerraLatY,bx,by);
% MerraLonX = MerraLonX(indexIN);
% MerraLatY = MerraLatY(indexIN);

mass = interp2(xp,yp,GRACEdiff,MerraLonX,MerraLatY)';

figure(10)
clf
scatter(totalmelt(indexIN'),-mass(indexIN'),5,MerraLonX(indexIN'),'x')
colorbar

figure(11)
clf
imagesc(MerraLonX(indexIN'))
colorbar
% temp = totalmelt.*indexIN';