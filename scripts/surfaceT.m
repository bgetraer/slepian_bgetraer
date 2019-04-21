%**************************************************************************
% Script for analyzing surface temperature over the Greenland Ice Sheet
% Merra2 Data from:
%
%   SCRIPT 1
%   Loads mean daily T values from netCDF files, isolates an estimate of
%   surface temperature across pressure levels, generates "melt day" data
%   defined as mean daily T < 0 degrees Celsius. Melt day anomalies per
%   month per year are calculated by removing the monthly mean for the
%   length of the record. Melt events over a duration (ie 3 consecutive
%   days) are also considered.
%
%
% Last modified by: bgetraer@princeton.edu 3/31/2019
%**************************************************************************

% locate slepian_bgetraer function and datafile directories, and set workspace
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
merra_dir = fullfile(dir,'datafiles/MERRA2');
addpath(fullfile(dir,'scripts/atmospheric_analysis/functions'))
addpath(fullfile(dir,'functions'))
setworkspace('/Users/benjamingetraer/Documents/IndependentWork/SH_Workspace');
datadir = fullfile(dir,'datafiles/');



ncdfDir = fullfile(merra_dir,'GrnlandSurfT');
matDir = fullfile(merra_dir,'MerraMat');
startdate = datenum(2003,01,01);
enddate = datenum(2017,12,31);
alldates = startdate:enddate;
%% LOAD NETCDF FILES AND SAVE THE SURFACE TEMPERATURE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Surface Skin Temperature for every day between 2003 and 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ST = load(fullfile(matDir,'TS2003-2017'),'data','t','spacelim');
ST.data = squeeze(ST.data);

T2M = load(fullfile(matDir,'T2M2003-2017'),'data','t','spacelim');
T2M.data = squeeze(T2M.data);

QV2M = load(fullfile(matDir,'QV2M2003-2017'),'data','t','spacelim');
QV2M.data = squeeze(QV2M.data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The GRIS Merra nan mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(datadir,'comparemassMerra'),'indexICEmerra','subREG');
GRISnan = double(indexICEmerra);
GRISnan(~indexICEmerra)=nan;


%% GET THE DATA ONLY OVER REGIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The GRIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRIS = struct;
GRIS.T2M = T2M.data.*GRISnan';
GRIS.ST = ST.data.*GRISnan';
GRIS.QV2M = QV2M.data.*GRISnan';

GRIS.T2Mvec = GRIS.T2M(~isnan(GRIS.T2M));
GRIS.STvec = GRIS.ST(~isnan(GRIS.ST));
GRIS.QV2Mvec = GRIS.QV2M(~isnan(GRIS.QV2M));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBREGIONS GREENLAND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REG = struct;
for i = 1:4
    REG.T2M{i} = T2M.data.*GRISnan'.*subREG.MERRA.nanindex{i}';
    REG.ST{i} = ST.data.*GRISnan'.*subREG.MERRA.nanindex{i}';
    REG.QV2M{i} = QV2M.data.*GRISnan'.*subREG.MERRA.nanindex{i}';
    
    REG.T2Mvec{i} = REG.T2M{i}(~isnan(REG.T2M{i}));
    REG.STvec{i} = REG.ST{i}(~isnan(REG.ST{i}));
    REG.QV2Mvec{i} = REG.QV2M{i}(~isnan(REG.QV2M{i}));
end
%% COMPARISON OF SKIN TEMP TO NEAR SURFACE TEMP ALL GREENLAND
nbins = 100;
%put the numbers into bins
[N,XEDGESgl,YEDGESgl,BINX,BINY] = histcounts2(GRIS.T2Mvec,GRIS.STvec,...
    nbins);
bin_meansgl = [];
for n = 1:nbins
    for m = 1:nbins
        bin_meansgl(n,m) = nanmean(GRIS.QV2Mvec(BINX == n & BINY == m))'; %average Y in each bin/subinterval
    end
end
bin_meansgl(isnan(bin_meansgl)) = 0;

%%
figure(1)
clf
hold on

histogram2('XBinEdges',XEDGESgl,'YBinEdges',YEDGESgl,'BinCounts',bin_meansgl,...
    'FaceColor','flat')
view(0,90)
% histogram2(GRIS.T2Mvec(:), GRIS.STvec(:),'FaceColor','flat')
% [a, ~, ~, ~, e]=pca([GRIS.T2Mvec(:) GRIS.STvec(:)]);
% mpca = a(1,1)/a(2,1);
% bpca = mean(GRIS.STvec(:)-mpca*GRIS.T2Mvec(:));
% pcline = plot3(xlim,mpca*xlim+bpca,[max(zlim)  max(zlim)]);
% set(pcline,'color','b','linewidth',2)

[m,f] = linear_m(GRIS.T2Mvec,GRIS.STvec);
slopeline = plot3(xlim,m(2).*xlim+m(1),[max(zlim)  max(zlim)]);
set(slopeline,'color','k','linewidth',2,'linestyle','--')

oneline = plot3(xlim,xlim,[max(zlim)  max(zlim)]);
set(oneline,'color','k','linewidth',2)

setST0 = plot3(xlim,[273.16 273.16],[max(zlim)  max(zlim)]);
set(setST0,'color','r','linewidth',1)
setT2M0 = plot3([273.16 273.16],ylim,[max(zlim)  max(zlim)]);
set(setT2M0,'color','r','linewidth',1)

ylabel('Surface Skin Temperature [K]')
xlabel('Near Surface 2M Air Temperature [K]')
% caxis([0 5E-3])
cb = colorbar;
ylabel(cb,'Near Surface 2M Specific Humidity [kg per kg]')
title('Greenland Ice Sheet')
axis square
legend([oneline slopeline],'1:1 line',...
    sprintf('m = %0.2f b = %0.2f',m(2),m(1)),...
    'location','northwest')


%% SUBREGIONS OF GREENLAND
% BIN THEM SO THAT YOU CAN SEE WHAT IS GOING ON
% THIS TAKES A WHILE

nbins = 100;

for i = 1:4
    fprintf('REGION %0.0d',i);
    %put the numbers into bins
    [N,XEDGES{i},YEDGES{i},BINX,BINY] = histcounts2(REG.T2Mvec{i},REG.STvec{i},...
        nbins);
    bin_means{i} = [];
    for n = 1:nbins
        for m = 1:nbins
            bin_means{i}(n,m) = nanmean(REG.QV2Mvec{i}(BINX == n & BINY == m))'; %average Y in each bin/subinterval
        end
    end
    bin_means{i}(isnan(bin_means{i})) = 0;
end
%%
ttltext = {'NORTHWEST','NORTHEAST','SOUTHEAST','SOUTHWEST'};
figure(3)
clf
for i = 1:4
    subplot(1,17,i*4-3:i*4)
    histogram2('XBinEdges',XEDGES{i},'YBinEdges',YEDGES{i},'BinCounts',bin_means{i},...
        'FaceColor','flat')
    view(0,90)
    ylim([200,290])
    xlim([200,290])
    [m,f] = linear_m(REG.T2Mvec{i},REG.STvec{i});
    hold on
    slopeline = plot3(xlim,m(2).*xlim+m(1),[1 1]);
    set(slopeline,'color','k','linewidth',2,'linestyle','--')
    
    oneline = plot3(xlim,xlim,[1 1]);
    set(oneline,'color','k','linewidth',2)
    
    setST0 = plot3(xlim,[273.16 273.16],[1 1]);
    set(setST0,'color','r','linewidth',1)
    setT2M0 = plot3([273.16 273.16],ylim,[1 1]);
    set(setT2M0,'color','r','linewidth',1)
    
    

    if i == 1
        ylabel('Surface Skin Temperature [K]')
%         set(gca,'yticklabels',yticks-273.16)
    else
        yticklabels ''
    end
    xlabel('Near Surface 2M Air Temperature [K]')
%     set(gca,'xticklabels',xticks-273.16)
    caxis([0 5E-3])
    title(ttltext{i})
    
    legend([oneline slopeline],'1:1 line',...
        sprintf('m = %0.2f b = %0.2f',m(2),m(1)),...
        'location','northwest')

    axis square
end
subplot(1,17,17)
xlim([200,290])
ylim([200,290])
axis off
caxis([0 5E-3])
cb = colorbar;
ylabel(cb,'Near Surface 2M Specific Humidity [kg per kg]')

%% PLOTTING DETAILS
oneline = refline(1,0);
set(oneline,'color','k','linewidth',2)

slopeline = refline(m(2),m(1));
set(slopeline,'color','k','linewidth',2,'linestyle','--')

setST0 = refline(0,273.16);
set(setST0,'color','r','linewidth',1)
setT2M0 = line([273.16 273.16],ylim);
set(setT2M0,'color','r','linewidth',1)

ylabel('Surface Skin Temperature [K]')
xlabel('Near Surface 2M Air Temperature [K]')
cb = colorbar;
ylabel(cb,'Near Surface 2M Specific Humidity [kg per kg]')

legend([slopeline],...
    sprintf('Least-squares linear regression: m = %0.2f b = %0.2f',m(2),m(1)),...
    'location','northwest')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DIVIDE BY SEASON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
season = alldates;
season(month(alldates)>=12 | month(alldates)<=2) = 1;
season(month(alldates)>=3 & month(alldates)<=5) = 2;
season(month(alldates)>=6 & month(alldates)<=8) = 3;
season(month(alldates)>=9 & month(alldates)<=11) = 4;

%% SEASONAL AIR TEMP ANOMALY

theseason = 4;


figure(4)
clf
col = [0.2422    0.1504    0.6603;...
    0.1540    0.5902    0.9218;...
    0.5044    0.7993    0.3480;...
    0.9769    0.0805   0.305];
sym = {'-^','-v','-<','->'};
seas = {'Winter (DJF)','Spring (MAM)','Summer (JJA)','Fall (SON)'};
for j = 1:4
    [t2m] = ...
        plotoverGRIS(alldates,REG.T2M{j});

    T2MA=anomMonth(t2m,alldates);
    
    [T2Mm,yr] = downsample_ts(T2MA(season==theseason),alldates(season==theseason),...
        'yearly');
    T2Mstd = downsample_ts(T2MA(season==theseason),alldates(season==theseason),'downsamplingperiod','yearly','function','nanstd');
   
    
    hold on
    h{j} = errorbar(yr,T2Mm,T2Mstd,'color',col(j,:),'linewidth',1.5);
    h{j} = plot(yr,T2Mm,sym{j},'color',col(j,:),'linewidth',1.5,...
        'markerfacecolor',col(j,:));

end

[t2m] = ...
    plotoverGRIS(alldates,GRIS.T2M);

T2MA=anomMonth(t2m,alldates);

[T2Mm,yr] = downsample_ts(T2MA,alldates,'yearly');
[m,f] = linear_m(yr,T2Mm);
h_line = refline(m(2),m(1));
set(h_line,'color','k','linewidth',1.5)


grid on
title(sprintf('%s 2M Air Temperature Anomaly',seas{theseason}))
ylabel('\circK')
lgd = legend([h{1},h{2},h{3},h{4},h_line],'nw','ne','se','sw',...
    sprintf('slope = %0.2f',m(2).*365.*10),'location','eastoutside');
title(lgd,{'mean 2MT','errorbars mark 1\sigma'})

xlim([datenum(2003,1,1),datenum(2017,12,31)])
    datetick('keeplimits')
%% EXAMINE AVERAGE OVER ICESHEET
figure(2)
[t2m,meanT2M,T2MA,fp,fl,ml,x,y] = plotoverGRIS(alldates,squeeze(T2M.data));

clf

[ts,meanTskin,TSA,fp,fl,ml,x,y] = plotoverGRIS(alldates,ST.data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,1)
hold on
ylabel('\circK')
title('mean surface skin T over GRIS')

subplot(3,3,4:5)
hold on
ylabel('\circK')
title('mean surface skin T anomaly over GRIS')

subplot(3,3,6)
ylabel('\circK')


ylabel('\circK')


%% LOAD PIXEL BY PIXEL GREENLAND TEMPERATURE OVER TIME
% [ Tslopes ] = linearMap( squeeze(T2M.data), alldates );

yrs = repmat((1980:2017),12,1);
months = repmat((1:12)',size(unique(yrs),1),1);
day = repmat(15,size(months,1),1);
allMdates = datenum(yrs(:),months,day);
imjset = squeeze(Temp.T2MMEAN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DIVIDE BY SEASON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
season = allMdates;
season(month(allMdates)>=12 | month(allMdates)<=2) = 1;
season(month(allMdates)>=3 & month(allMdates)<=5) = 2;
season(month(allMdates)>=6 & month(allMdates)<=8) = 3;
season(month(allMdates)>=9 & month(allMdates)<=11) = 4;

% [ TslopesFull ] = linearMap( imjset, allMdates );
% [ TslopesGRACE ] = linearMap( imjset(:,:,year(allMdates)>=2003),...
%     allMdates(year(allMdates)>=2003));
% [ TslopesGRACE10,~,R ] = linearMap( imjset(:,:,year(allMdates)>=2003& year(allMdates)<=2012),...
%     allMdates(year(allMdates)>=2003& year(allMdates)<=2012));
% [ TslopesGRACE5 ] = linearMap( imjset(:,:,year(allMdates)>=2013& year(allMdates)<=2017),...
%     allMdates(year(allMdates)>=2013& year(allMdates)<=2017));

% ONLY SUMMER/WINTER

[ TslopesFullSummer,~, Rsummer ] = linearMap( imjset(:,:,season==3), allMdates(season==3) );
% [ TslopesGRACESummer ] = linearMap( imjset(:,:,season==3 & year(allMdates)>=2003),...
%     allMdates(season==3 & year(allMdates)>=2003));
% [ TslopesGRACESummer10, ~, R] = linearMap( imjset(:,:,season==3 & year(allMdates)>=2003& year(allMdates)<=2012),...
%     allMdates(season==3 & year(allMdates)>=2003& year(allMdates)<=2012));
% [ TslopesGRACESummer5 ] = linearMap( imjset(:,:,season==3 & year(allMdates)>=2013& year(allMdates)<=2017),...
%     allMdates(season==3 & year(allMdates)>=2013& year(allMdates)<=2017));

% [ TslopesFullWinter ] = linearMap( imjset(:,:,season==1), allMdates(season==1) );
% [ TslopesGRACEWinter ] = linearMap( imjset(:,:,season==1 & year(allMdates)>=2003),...
%     allMdates(season==1 & year(allMdates)>=2003));
% [ TslopesGRACEWinter10 ] = linearMap( imjset(:,:,season==1 & year(allMdates)>=2003& year(allMdates)<=2012),...
%     allMdates(season==1 & year(allMdates)>=2003& year(allMdates)<=2012));

% [ TslopesFullnotSummer ] = linearMap( imjset(:,:,season~=3), allMdates(season~=3) );
% [ TslopesGRACEnotSummer ] = linearMap( imjset(:,:,season~=3 & year(allMdates)>=2003),...
%     allMdates(season~=3 & year(allMdates)>=2003));
% [ TslopesGRACEnotSummer10 ] = linearMap( imjset(:,:,season~=3 & year(allMdates)>=2003& year(allMdates)<=2012),...
%     allMdates(season~=3 & year(allMdates)>=2003& year(allMdates)<=2012));
% [ TslopesGRACEnotSummer5 ] = linearMap( imjset(:,:,season~=3 & year(allMdates)>=2013& year(allMdates)<=2017),...
%     allMdates(season~=3 & year(allMdates)>=2013& year(allMdates)<=2017));

% save(fullfile(matDir,'tslopes'),'TslopesFull','TslopesFullnotSummer',...
%     'TslopesFullSummer','TslopesGRACE10','TslopesGRACEnotSummer10',...
%     'TslopesGRACESummer10','TslopesGRACE','TslopesGRACEnotSummer',...
%     'TslopesGRACESummer','TslopesGRACESummer5');

%% SLOPE IMAGE
ax5 = figure(5);
clf
sth = suptitle('T2M Trends Over The Greenland Ice Sheet [\circC per decade]');
sth.Position = sth.Position + [0 0.04 0];
set(gcf,'Position',[1 91 1089 607])
imj = {TslopesFull,TslopesFullnotSummer,TslopesFullSummer,...
    TslopesGRACE10,TslopesGRACEnotSummer10,TslopesGRACESummer10,...
    TslopesGRACE,TslopesGRACEnotSummer,TslopesGRACESummer10};
capY = {'1980-2017','2003-2012','2003-2017'};
capS = {'All Seasons','Fall/Winter/Spring','Summer'};
for i = 1:9
    % PLOT
    if i < 4
        subplot(3,4,i)
    elseif i < 7 
        subplot(3,4,i+1)
    else
        subplot(3,4,i+2)
    end
    hold on
    thisimj = interp2(pM.X,pM.Y,imj{i}'.*365.*10,GRACEX,GRACEY);
    
    fill(pG.gx./resamprate+0.5,pG.gy./resamprate+0.5,[0.6 0.6 0.6],'EdgeColor','none')
    hold on
    plot(pG.bx./resamprate+0.5,pG.by./resamprate+0.5,'--k','linewidth',0.2)
    set(gca,'ydir','reverse')
    
    imagesc(thisimj,'AlphaData',~isnan(thisimj));
%     if any(i == 1:2)
%     end
    axis off
    colormap(ax5,jet(20))
    caxis([prctile(thisimj(:),1), prctile(thisimj(:),99)])
    cb = colorbar;
    if abs(diff(minmax(cb.Ticks)))> 1
        cb.Ticks = [-2:0.4:2];
    elseif abs(diff(minmax(cb.Ticks)))>0.5
        cb.Ticks = [-2:0.2:2];
    else
        cb.Ticks = [-2:0.05:2];
    end
    cb.Position = cb.Position + 1e-10;
%     ylabel(cb, '\circC per decade');
%     if any(i == 1:3)
%         capy = capY{1};
%     elseif any(i == 4:6)
%         capy = capY{2};
%     else
%         capy = capY{3};
%     end
%     if any(i == [1 4 7])
%         caps = capS{1};
%     elseif any(i == [2 5 8])
%         caps = capS{2};
%     else
%         caps = capS{3};
%     end
%     th = title(sprintf('T2M Trend %s, %s',capy,caps))
    if any(i == 1:3)
        ttlh = title(capS{i});
        ttlh.Position = ttlh.Position+[0 -10 0];
    end
    
    if any(i == [1 4 7])
        txt = sprintf('\\bf{%s}',capY{i==[1 4 7]});
        txth = text(0,100,txt,'fontsize',11);
        txth.Rotation = 90;
    end
end


subplot(3,4,4)
hold on
scatter(TslopesFull(:),TslopesFullSummer(:),'x')
scatter(TslopesFull(:),TslopesFullnotSummer(:),'+')

% scatter(Aall,Asum,'x')
% scatter(Aall,Anot,'+')
axis([-0.73 1 -0.73 1]/3650)
axis square %off
refline(1,0)
set(gca,'xticklabels',round(xticks*365*10,2),'yticklabels',round(yticks*365*10,2))
xlabel('all seasons')
ylabel('summer/non-summer')

subplot(3,4,8)
hold on
scatter(TslopesGRACE10(:),TslopesGRACESummer10(:),'x')
scatter(TslopesGRACE10(:),TslopesGRACEnotSummer10(:),'+')
axis([-4 2.2 -4 2.2]/3650)
axis square %off
refline(1,0)
set(gca,'xticklabels',xticks*365*10,'yticklabels',round(yticks*365*10,2))

xlabel('all seasons')
ylabel('summer/non-summer')

subplot(3,4,12)
hold on
scatter(TslopesGRACE(:),TslopesGRACESummer(:),'x')
scatter(TslopesGRACE(:),TslopesGRACEnotSummer(:),'+')
axis([-2.2 0.8 -2.2 0.8]/3650)
axis square %off
refline(1,0)
refline(1,0)
set(gca,'xticklabels',xticks*365*10,'yticklabels',round(yticks*365*10,2))

xlabel('all seasons')
ylabel('summer/non-summer')

lgd = legend('all vs summer','all vs non-summer');
% lgd.Orientation = 'horizontal';
lgd.Position = [0.8432 0.1415 0.0902 0.0896];
%%
figure(6)
clf
hold on
% 
% Gsum = (TslopesGRACESummer(:)-nanmean(TslopesGRACESummer(:)))./nanstd(TslopesGRACESummer(:));
% Gall = (TslopesGRACE(:)-nanmean(TslopesGRACE(:)))./nanstd(TslopesGRACE(:));
% Gnot = (TslopesGRACEnotSummer(:)-nanmean(TslopesGRACEnotSummer(:)))./nanstd(TslopesGRACEnotSummer(:));
% 
% Asum = (TslopesFullSummer(:)-nanmean(TslopesFullSummer(:)))./nanstd(TslopesFullSummer(:));
% Aall = (TslopesFull(:)-nanmean(TslopesFull(:)))./nanstd(TslopesFull(:));
% Anot = (TslopesFullnotSummer(:)-nanmean(TslopesFullnotSummer(:)))./nanstd(TslopesFullnotSummer(:));
% 
% scatter(Asum,Aall)
% scatter(Asum,Gsum)
% scatter(Anot,Aall)
% scatter(Gnot,Gall)
% scatter(TslopesFullnotSummer(:),TslopesFull(:))
scatter(TslopesGRACESummer(:),TslopesFullSummer(:))
% scatter(TslopesGRACESummer(:),TslopesGRACE(:))
% axis([-1 2.5 -1 2.5]*1E-4)
axis square
refline(1,0)

%% PLOT MEAN T PER REGION
load(fullfile(datadir,'comparemassMerra'),'subREG','indexICE','indexICEmerra','indexICEnan',...
    'indexGL','indexGLnan','GRACEX','GRACEY','resamp','resamprate')

figure(6)
clf
hold on
for i = 1:4
    nanindex = double(subREG.MERRA.index{i}');
    nanindex(~nanindex) = nan;
    thisregvec = reshape(ST.data.*nanindex,...
        size(ST.data,1)*size(ST.data,2),size(ST.data,3));
    thismeltvec = thisregvec>273.15;
    thisregmean = nanmean(thisregvec,1);
    plot(alldates,thisregmean)
end
grid on
%% MELT DAYS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Classify melt day for every day between 1987 and 2018:
%   meanT(day) >    273.16 K    --> 1
%   meanT(day) <=   273.16 K    --> 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meltThreshold = 273.16;
meltThreshold = 270;
allmeltMap = T2M.data >= meltThreshold;

%% SHOW TOTAL MELT DAYS OVER GRACE
GRACEmeltMap = allmeltMap(:,:,find(firstdate:enddate == startdate):end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOTAL MELT DAYS INSIDE OF GREENLAND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalmelt = sum(allmeltMap,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THE COORDINATE SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(matDir,'projectMERRAGL'),'Flon2x','Flat2y','gx','gy','bx',...
    'by','X','Y','LON','LAT','contx','conty','azores','reykjavik');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT TOTAL MELT DAYS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
clf
GRISmerranan = double(indexICEmerra);
GRISmerranan(~indexICEmerra) = nan;
imj = totalmelt./size(allmeltMap,3).*GRISmerranan';
imPlot(imj)
hold on
plot(bx,by,'k--')
plot(gx,gy,'k')
title('FRACTION OF DAYS WITH MEAN(T) > -4.16 C\circ JAN2003-JAN2018')
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
    imj = avgMonthMelt(:,:,m).*GRISmerranan';
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
    
    subplot(3,5,i)
    
    imj = thisYMelt;
    
    thisimj = interp2(pM.X,pM.Y,imj',GRACEX,GRACEY);
    
    fill(pG.gx./resamprate+0.5,pG.gy./resamprate+0.5,[0.6 0.6 0.6],'EdgeColor','none')
    hold on
    plot(pG.bx./resamprate+0.5,pG.by./resamprate+0.5,':k','linewidth',0.2)
    set(gca,'ydir','reverse')
    
    imagesc(thisimj,'AlphaData',~isnan(GRISnan));
    
%     imPlot(imj,parula,cax)
    hold on
%     plot(bx,by,'k--')
%     plot(gx,gy,'k')
    title(num2str(thisYear))
    caxis(cax)
%     colorbar('eastoutside')
    axis tight  off
end

suptitle('Melt days (daily mean T2M > 270^{\circ}K) anomaly by year, normalized by monthly mean')

axes('Position',[0.9 0.3 0.02 0.4])
axis off
caxis(cax)
cb = colorbar;
cb.Position = cb.Position + [0 0 0.01 0];
ylabel(cb,'# days, normalized by monthly mean','fontsize',10)
