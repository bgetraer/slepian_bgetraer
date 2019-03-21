%**************************************************************************
% Script for analyzing temperature over the Greenland Ice Sheet
% Merra2 Data from:
% https://disc.gsfc.nasa.gov/datasets/M2I6NPANA_V5.12.4/summary?keywords=%22MERRA-2%22
%
%   SCRIPT 1
%   Loads mean daily T values from netCDF files, isolates an estimate of
%   surface temperature across pressure levels, generates "melt day" data
%   defined as mean daily T < 0 degrees Celsius. Melt day anomalies per
%   month per year are calculated by removing the monthly mean for the
%   length of the record. Melt events over a duration (ie 3 consecutive
%   days) are also considered.
%
%   NEXT: COMPARETEMPMASS
%   
% SEE ALSO:
%   GETMERRA2VAR, FLOORDATA, SURFACETBYYEAR, MELTMAP, MELTTHRESH, 
%   PROJECTMERRA, ANOMMONTH
%
% Last modified by: bgetraer@princeton.edu 3/13/2019
%**************************************************************************

% locate slepian_bgetraer function and datafile directories, and set workspace
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
merra_dir = fullfile(dir,'datafiles/MERRA2');
addpath(fullfile(dir,'scripts/atmospheric_analysis/functions'))
addpath(fullfile(dir,'functions'))
setworkspace('/Users/benjamingetraer/Documents/IndependentWork/SH_Workspace');


%% LOAD NETCDF FILES AND SAVE THE "TEMPERATURE FLOOR"
ncdfDir = fullfile(merra_dir,'MerraGrnlandTMEAN');
matDir = fullfile(merra_dir,'MerraMat');
firstdate = datenum(2003,01,01);
startdate = datenum(2003,01,01);
enddate = datenum(2018,01,01);
alldates = startdate:enddate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%n%%%%%%%%%%
% Get the Surface Temperature for every day between 1987 and 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[surfTFile, surfTyear, thespacelim] = surfaceTbyYear(firstdate, enddate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Classify melt day for every day between 1987 and 2018:
%   meanT(day) >    272.15 K    --> 1
%   meanT(day) <=   272.15 K    --> 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[meltmapFile] = meltMap(firstdate, enddate);

%% SHOW TOTAL MELT DAYS OVER GRACE
load(meltmapFile)
GRACEmeltMap = allmeltMap(:,:,find(firstdate:enddate == startdate):end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOTAL MELT DAYS INSIDE OF GREENLAND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalmelt = sum(allmeltMap,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THE COORDINATE SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ Flon2x,Flat2y,gx,gy,bx,by,contx,conty,azores,reykjavik, X,Y,LON,LAT] = ...
    projectMERRA( totalmelt, thespacelim );

save(fullfile(matDir,'projectMERRAGL'),'Flon2x','Flat2y','gx','gy','bx',...
    'by','X','Y','LON','LAT','contx','conty','azores','reykjavik');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MASK OF GREENLAND IN IMAGE BASIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz = size(X);
maskarray = inpolygon(X(:),Y(:),bx,by);         % array of points inside Greenland
mask = imfill(1*reshape(maskarray,sz)');           % matrix of points inside Greenland
mask(mask==0) = nan;
totalmelt = totalmelt.*mask;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT TOTAL MELT DAYS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
clf
imj = totalmelt./size(allmeltMap,3);
imPlot(imj,'jet')
hold on
plot(bx,by,'k--')
plot(gx,gy,'k')
title('FRACTION OF DAYS WITH MEAN(T) > 0 C\circ JAN2003-JAN2018')
latlonaxis(xlim,ylim,thespacelim,10)
colorbar('eastoutside')
axis tight

%% CALCULATE AND PLOT MONTHLY ANOMALY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CALCULATE MONTHLY ANOMALY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ normMeltMap, avgMonthMelt ] = anomMonth( allmeltMap, alldates);

save(fullfile(matDir,'normalizedMeltMap'),'normMeltMap','avgMonthMelt',...
    'alldates','allmeltMap');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT MONTHLY ANOMALY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
clf
cax = [0 1];

for m = 1:12
    subplot(4,3,m)
    imj = avgMonthMelt(:,:,m).*mask;
    imPlot(imj,'jet',cax)
    hold on
    plot(bx,by,'k--')
    plot(gx,gy,'k')
    title(datestr(datenum(1,m,1),'mmmm'))
    colorbar('eastoutside')
    axis tight
    axis off
end

y = unique(year(alldates));

suptitle(sprintf('Average fraction of melt days, %i-%i',...
    min(y),max(y)))


%% MAP OUT ALL YEARS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT MELT DAYS ANOMALY BY YEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
clf
cax = [-25 30];

for i = 1:length(y)
    thisYear = y(i);
    thisYMelt = sum(normMeltMap(:,:,year(alldates)==thisYear),3);

    total(i) = sum(thisYMelt(:));
    subplot(4,4,i)
    imj = thisYMelt.*mask;
    imPlot(imj,parula,cax)
    hold on
    plot(bx,by,'k--')
    plot(gx,gy,'k')
    title(num2str(thisYear))
    colorbar('eastoutside')
    axis tight
    axis off
end

suptitle('Melt days anomaly by year, normalized by monthly mean')

%% MELT EVENTS OVER A DURATION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CALCULATE 3 DAY MELT EVENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
past3days = movsum(allmeltMap,[2 0],3);
past3days(past3days<3) = 0;
past3days(past3days>=3) = past3days(past3days>=3)-2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NORMALIZE BY MONTHLY MEAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ normMelt3, avgMonthMelt3 ] = anomMonth( past3days, alldates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT MONTHLY ANOMALY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
clf
cax = [0 1];

for m = 1:12
    subplot(4,3,m)
    imj = avgMonthMelt3(:,:,m).*mask;
    imPlot(imj,'jet',cax)
    hold on
    plot(bx,by,'k--')
    plot(gx,gy,'k')
    title(datestr(datenum(1,m,1),'mmmm'))
    colorbar('eastoutside')
    axis tight
    axis off
end

y = unique(year(alldates));

suptitle(sprintf('Average fraction of 3 consecutive melt days, %i-%i',...
    min(y),max(y)))

%% MAP OUT ALL YEARS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT 3 CONSECUTIVE MELT DAYS ANOMALY BY YEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
clf
cax = [-25 30];

for i = 1:length(y)
    thisYear = y(i);
    thisYMelt = sum(normMelt3(:,:,year(alldates)==thisYear),3);

    subplot(4,4,i)
    imj = thisYMelt.*mask;
    imPlot(imj,'parula',cax)
    hold on
    plot(bx,by,'k--')
    plot(gx,gy,'k')
    title(num2str(thisYear))
    colorbar('eastoutside')
    axis tight
    axis off
end

suptitle('3 Consecutive melt days anomaly by year, normalized by monthly mean')



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



datadir = fullfile(dir,'datafiles');
addpath(dir,datadir);

load ptsGL 
load im_tools 
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
%%


%%


GRACEdiff = D(:,:,monthnum(6,2006,thedates)) - D(:,:,1);
wavename = 'haar';

level = wmaxlev(size(GRACEdiff),wavename);
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

% plot(MerraLonX,MerraLatY,'.')

figure(2)
clf
imagesc(GRACEdiffprime,'AlphaData',~isnan(GRACEdiffprime));
hold on
axis image
plot(gx,gy,'k-')
colormap(bluewhitered([],1))
axis off


%%

[ Flon2x,Flat2y,gx,gy,bx,by,contx,conty ] = projectMERRA( totalmelt, thespacelim );
figure(1)
clf
    imPlot(totalmelt)
    hold on
    plot(bx,by,'w:')
    plot(contx,conty,'k-','linewidth',1)


%%
%   make it into a function that looks for the file then creates it if not
%   there

%   create a function that says if a certain point is melting in a given
%   date


figure(1)
clf
for i = 1:size(meltmap,3)
    imPlot(meltmap(:,:,i))
    hold on
    plot(bx,by,'w:')
    plot(contx,conty,'k-','linewidth',1)
    i
    pause(0.01)
end
%%

[ Flon2x,Flat2y,gx,gy,bx,by,contx,conty ] = projectMERRA( tdata, thespacelim );
intMeltdays = zeros(size(tdata,1),size(tdata,2));

figure(1)
for i = 1:size(tdata,4)
    t_floor = floorData(tdata(:,:,:,i));    
    
    clf
    subplot(1,3,1)
    imPlot(t_floor);
        hold on

    plot(bx,by,'w:')
    plot(contx,conty,'k-','linewidth',1)
    colorbar
    latlonaxis( xlim, ylim, thespacelim, 5 );
    
    % put the date on
    text(Flon2x(-75),Flat2y(85),thetimedata(i),'color','k',...
        'verticalalignment','top');
    
    subplot(1,3,2)
    
    % plotting
    meltmap = isMelting(t_floor);
    imPlot(meltmap,bone);
    hold on
    plot(bx,by,'w:')
    plot(contx,conty,'k-','linewidth',1)

%     colorbar
    %axis
    latlonaxis( xlim, ylim, thespacelim, 5 );
    
    
    % put the date on
    text(Flon2x(-89),Flat2y(89),thetimedata(1),'color','w');
    
    pause(0.1)
    
    subplot(1,3,3)
    % plotting
    
    intMeltdays = intMeltdays+meltmap;
    
    imPlot(intMeltdays,bone);
    hold on
    plot(bx,by,'w:')
    plot(contx,conty,'k-','linewidth',1)

%     colorbar
    %axis
    latlonaxis( xlim, ylim, thespacelim, 5 );
    
    % put the date on
    text(Flon2x(-89),Flat2y(89),thetimedata(1),'color','w');
    
%     pause(0.1)
    
end