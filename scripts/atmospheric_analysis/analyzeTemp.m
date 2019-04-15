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
datadir = fullfile(dir,'datafiles/');


%% LOAD NETCDF FILES AND SAVE THE "TEMPERATURE FLOOR"
ncdfDir = fullfile(merra_dir,'MerraGrnlandT');
matDir = fullfile(merra_dir,'MerraMat');
firstdate = datenum(2003,01,01);
startdate = datenum(2003,01,01);
enddate = datenum(2017,12,31);
alldates = startdate:enddate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%n%%%%%%%%%%
% Get the Surface Temperature for every day between 1987 and 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [surfTFile, surfTyear, thespacelim] = surfaceTbyYear(firstdate, enddate,...
%     ncdfDir, matDir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Classify melt day for every day between 1987 and 2018:
%   meanT(day) >    272.15 K    --> 1
%   meanT(day) <=   272.15 K    --> 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [meltmapFile] = meltMap(firstdate, enddate);
T2M = load(fullfile(matDir,'T2M2003-2017'),'data','t','spacelim');
T2M.data = squeeze(T2M.data);
%% SHOW TOTAL MELT DAYS OVER GRACE


allmeltMap = T2M.data >= 273.16;

GRACEmeltMap = allmeltMap(:,:,find(firstdate:enddate == startdate):end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOTAL MELT DAYS INSIDE OF GREENLAND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalmelt = sum(allmeltMap,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THE COORDINATE SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ Flon2x,Flat2y,gx,gy,bx,by,contx,conty,azores,reykjavik, X,Y,LON,LAT] = ...
    projectMERRA( totalmelt, T2M.spacelim );

save(fullfile(matDir,'projectMERRAGL'),'Flon2x','Flat2y','gx','gy','bx',...
    'by','X','Y','LON','LAT','contx','conty','azores','reykjavik');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MASK OF GREENLAND IN IMAGE BASIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz = size(X);
maskarray = inpolygon(X(:),Y(:),bx,by);         % array of points inside Greenland
mask = imfill(1*reshape(maskarray,sz)');           % matrix of points inside Greenland
mask(mask==0) = nan;
totalmelt = totalmelt.*mask;

% MERRA projection data
pM = load(fullfile(matDir,'projectMERRAGL'));
subregions = load(fullfile(datadir,'subregions'),'LONBUF','LATBUF',...
    'lonICE','latICE','AREA');
ice = struct;
ice.MERRA.X = pM.Flon2x(subregions.lonICE-360);
ice.MERRA.Y = pM.Flat2y(subregions.latICE);
indexICE = inpolygon(pM.X,pM.Y,ice.MERRA.X,ice.MERRA.Y);
indexICEnan = double(indexICE);
indexICEnan(indexICEnan==0) = nan;
%% SHOW TOTAL GREENLAND TEMPERATURE OVER TIME
yr = 2003:2017;
meanT = [];
for i = 1:length(yr)
    load(surfTFile{i})
    meanT = [meanT; squeeze(nanmean(nanmean(floorT.*GRISnan',1),2))];
    
    clear floorT
end
meanT = meanT(1:length(alldates));
meanTmean = mean(meanT);
meanT = meanT-meanTmean;
%% PLOT MEAN T OVER GREENLAND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DIVIDE BY SEASON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
season = alldates;
season(month(alldates)>=1 & month(alldates)<=2) = 1;
season(month(alldates)>=3 & month(alldates)<=5) = 2;
season(month(alldates)>=6 & month(alldates)<=9) = 3;
season(month(alldates)>=10 & month(alldates)<=12) = 4;
theseason = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DETREND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mp,fp,Rsq,power_m] = periodic_m(alldates,meanT,[365,365/2]);
meanTA = meanT-fp;

[ml,fl,DELTA,mint,t] = linear_m(alldates(season==theseason),meanTA(season==theseason));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
clf
subplot(3,1,1)
hold on
% plot(alldates(year(alldates)<2017), meanT(year(alldates)<2017)+meanTmean,'-')
plot(alldates, meanT+meanTmean,'-')
plot(alldates, fp+meanTmean,'-','linewidth',1)
datetick
grid on
axis tight
title('mean surface T over GRIS')

subplot(3,3,4:5)
hold on
cmap = repmat([0 0.4470 0.7410],4,1);
cmap(theseason,:) = 0;
colormap(cmap)
scatter(alldates,meanTA,1,season,'x')
% plot(alldates, meanTA,'-')
plot(alldates(season==theseason), fl,'-r','linewidth',1)
datetick
grid on
axis tight
title('mean surface T anomaly over GRIS')


subplot(3,3,6)
hold on
xoff = 0;
x = 1:365;
y = fp(x+xoff);
plot(x,y+meanTmean,'-','linewidth',1.5)
grid on
axis tight
datetick('x','mmm')
title('seasonal cycle removed')

subplot(3,1,3)
hold on

BP = zeros(120,length(yr))+nan;

seasonmean = mean(meanTA(season==theseason)+meanTmean);
for i = 1:length(yr)
    
    yeardata = meanTA(season==theseason & year(alldates)==yr(i))+meanTmean;
    BP(1:length(yeardata),i) = yeardata-seasonmean;
end

boxplot(BP,yr,'symbol','rx')
grid on
title('JJAS summer anomalies')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD PIXEL BY PIXEL GREENLAND TEMPERATURE OVER TIME
yr = 2003:2017;
iceT = T2M.data;

%% NORMALIZE PIXEL BY PIXEL
iceTvec = reshape(iceT,size(iceT,1)*size(iceT,2),size(iceT,3));
meaniceTvec = nanmean(iceTvec,2);

for i = 1:size(iceT,1)*size(iceT,2)
    if ~isnan(meaniceTvec(i))
        
        [mp,fp(i,:)] = periodic_m(alldates,iceTvec(i,:)-meaniceTvec(i),[365,365/2]);
        
        [ml(:,i)] = linear_m(alldates(season==theseason),iceTvec(i,(season==theseason))-meaniceTvec(i)-fp(i,(season==theseason)));
    else
        ml(:,i) = [0;0];
    end
end
%% 
n = 4879;
Tslopes = ml(2,:);
% Tslopes(n) = 1E6;
Tslopes = reshape(Tslopes,107,53);%.*365.*10;
% Tslopes = reshape(meaniceTvec,107,53);
figure(4)
clf
% imPlot(imgaussfilt(Tslopes,2).*indexICEnan')
% colorbar
plot(alldates,iceTvec(n,:))
% hold on 
plot(alldates,fp(n,:)+meaniceTvec(n))
y = iceTvec(n,:)-fp(n,:)-meaniceTvec(n);
plot(alldates,y)
[mlin] = linear_m(alldates,y);
axis tight 
datetick

figure(5)
clf
imj = imgaussfilt(Tslopes,2).*indexICEnan';
imPlot(imj)
colorbar

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
latlonaxis(xlim,ylim,T2M.spacelim,10)
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
