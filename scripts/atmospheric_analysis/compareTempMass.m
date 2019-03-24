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

%% Load data to structures
meltData = load(fullfile(matDir,'normalizedMeltMap'));
% MERRA projection data
pM = load(fullfile(matDir,'projectMERRAGL'));
% GRACE projection data
pG = load('im_tools');
ptsGL = load('ptsGL');
% GRACE mass data
massData = load('im_seqSH');
% subregion outlines
subregions = load(fullfile(datadir,'subregions'),'LONBUF','LATBUF');


%% PROJECT GRACE ONTO MERRA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resample the grid from the cubed sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resamp = 1:2:256;
resampleLON = ptsGL.lond(resamp,resamp);
resampleLAT = ptsGL.latd(resamp,resamp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transform coordinates from Lat/Lon to the MERRA X/Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRACEX = pM.Flon2x(resampleLON-360);
GRACEY = pM.Flat2y(resampleLAT);

subREG = struct;
subREG.XBUF = cell(1,4);
subREG.YBUF = cell(1,4);
for i = 1:length(subregions.LATBUF)
    subREG.XBUF{i} = pM.Flon2x(subregions.LONBUF{i}-360);
    subREG.YBUF{i} = pM.Flat2y(subregions.LATBUF{i});
end
summitLONLAT = [-37.58,72.57];
summit(1) = pM.Flon2x(-37.58);
summit(2) = pM.Flat2y(72.57);


%% Index by region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Index coordinated INSIDE of buffered greenland
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indexGL = inpolygon(GRACEX,GRACEY,pM.bx,pM.by);
indexGLmerra = inpolygon(pM.X,pM.Y,pM.bx,pM.by);

subREG.index  = cell(1,4);
for i = 1:length(subREG.index)
    subREG.index{i} = inpolygon(GRACEX,GRACEY,subREG.XBUF{i},subREG.YBUF{i});
end
%% COMPARE TOTALS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get total maps for GRACE and MERRA (filtered)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalmelt = imgaussfilt(sum(meltData.allmeltMap,3),0.1);
GRACEdiff = massData.D(:,:,end) - massData.D(:,:,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resample the maps for GRACE and MERRA (filtered)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mass = GRACEdiff(resamp,resamp);
mass = mass.*indexGL;

melt = interp2(pM.X,pM.Y,totalmelt',GRACEX,GRACEY);
melt = melt.*indexGL;

%% DISTANCE BETWEEN CENTER AND EDGE OF BUFFER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE A RELATIVE DISTANCE FROM EDGE OF ICESHEET TO SUMMIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance to summit
d2c = distance(resampleLAT,resampleLON,summitLONLAT(2),summitLONLAT(1));
% distance to edge of buffer
sz = size(GRACEX);
[~, d2b] = dsearchn([pM.bx,pM.by],[GRACEX(:),GRACEY(:)]);
d2b = reshape(d2b,sz);
% relative distance measure
distweight = d2b./(d2c + d2b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the relative distance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
clf
hold on
indNan = double(indexGL);
indNan(~indNan) = nan;
surf(GRACEX.*indNan,GRACEY.*indNan,distweight.*indNan);
plot3(summit(1),summit(2),2,'^r','markerfacecolor','r')
cmap = colormap(jet);
fill(pM.bx,pM.by,cmap(1,:))
colorbar
axis tight off
set(gca,'ydir','reverse')

%%
figure(11)
clf
    scatter(melt(:).*indexGL(:),-mass(:).*indexGL(:),5,distweight(:).*indexGL(:),'x')
colormap(jet)
colorbar

%% COMPARE TOTAL TEMP AND MASS FOR DIFFERENT SIDES OF GREENLAND

figure(11)
clf
suptitle('Mass loss to melt day comparison, 2003-2017')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAST AND WEST MASS LOSS v MELTDAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
region = {[1 4],[2 3]};
for j = 1:2
    subplot(2,2,j)
    hold on;grid on
    for reg = region{j}
        scatter(melt(:).*subREG.index{reg}(:),-mass(:).*subREG.index{reg}(:),5,distweight(:).*subREG.index{reg}(:),'x')
    end
    colormap(jet)
    xlabel('number of melt days')
    ylabel('Kg per m^2 of mass loss')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EAST AND WEST REGIONS AND DISTANCE COLORING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)
hold on
indNan = double(indexGL);
indNan(~indNan) = nan;
surf(GRACEX.*indNan,GRACEY.*indNan,distweight.*indNan);
plot3(summit(1),summit(2),2,'^k','markerfacecolor','r')
cmap = colormap(jet);
fill(pM.bx,pM.by,cmap(1,:))
axis tight off
set(gca,'ydir','reverse')
fill3(subREG.XBUF{1},subREG.YBUF{1},repmat(2,size(subREG.XBUF{1})),'w','FaceAlpha',0.5)
fill3(subREG.XBUF{4},subREG.YBUF{4},repmat(2,size(subREG.XBUF{4})),'w','FaceAlpha',0.5)


subplot(2,2,4)
hold on
surf(GRACEX.*indNan,GRACEY.*indNan,distweight.*indNan);
plot3(summit(1),summit(2),2,'^k','markerfacecolor','r')
cmap = colormap(jet);
fill(pM.bx,pM.by,cmap(1,:))
axis tight off
set(gca,'ydir','reverse')
fill3(subREG.XBUF{2},subREG.YBUF{2},repmat(2,size(subREG.XBUF{2})),'w','FaceAlpha',0.5)
fill3(subREG.XBUF{3},subREG.YBUF{3},repmat(2,size(subREG.XBUF{3})),'w','FaceAlpha',0.5)


%%

load passIndexbyArea
immask = load('GLimagemask');

figure(5)
clf
immask.A(~immask.A)=nan;
imagesc(GRACEdiff.*immask.A)
colormap(bluewhitered)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROJECT MERRA ONTO GRACE
% Get the MERRA LAT LON grid 
MerraLonX = pG.Fx(pM.LON,pM.LAT);
MerraLatY = pG.Fy(pM.LON,pM.LAT);
indexGL = inpolygon(MerraLonX,MerraLatY,pG.bx,pG.by);
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
