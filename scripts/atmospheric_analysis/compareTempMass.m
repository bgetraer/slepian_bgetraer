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
subregions = load(fullfile(datadir,'subregions'),'LONBUF','LATBUF',...
    'lonICE','latICE','AREA');

%% PROJECT GRACE ONTO MERRA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resample the grid from the cubed sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resamprate = 2;
resamp = 1:resamprate:256;
resampleLON = ptsGL.lond(resamp,resamp);
resampleLAT = ptsGL.latd(resamp,resamp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transform coordinates from Lat/Lon to the MERRA X/Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRACEX = pM.Flon2x(resampleLON-360);
GRACEY = pM.Flat2y(resampleLAT);

subREG = struct;
% subREG.MERRA.XBUF = cell(1,4);
% subREG.MERRA.YBUF = cell(1,4);
subREG.GRACE.XBUF = cell(1,4);
subREG.GRACE.YBUF = cell(1,4);
for i = 1:length(subregions.LATBUF)
    subREG.MERRA.XBUF{i} = pM.Flon2x(subregions.LONBUF{i}-360);
    subREG.MERRA.YBUF{i} = pM.Flat2y(subregions.LATBUF{i});
    subREG.GRACE.XBUF{i} = pG.Fx(subregions.LONBUF{i},subregions.LATBUF{i});
    subREG.GRACE.YBUF{i} = pG.Fy(subregions.LONBUF{i},subregions.LATBUF{i});
end

ice = struct;
ice.MERRA.X = pM.Flon2x(subregions.lonICE-360);
ice.MERRA.Y = pM.Flat2y(subregions.latICE);
ice.GRACE.X = pG.Fx(subregions.lonICE,subregions.latICE);
ice.GRACE.Y = pG.Fy(subregions.lonICE,subregions.latICE);

summitLONLAT = [-37.58+360,72.57];
summit(1) = pG.Fx(summitLONLAT(1),summitLONLAT(2));
summit(2) = pG.Fy(summitLONLAT(1),summitLONLAT(2));

%% Index by region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Index coordinates inside of buffered greenland
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indexGL = inpolygon(GRACEX,GRACEY,pM.bx,pM.by);
indexGLnan = double(indexGL);
indexGLnan(indexGLnan==0) = nan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Index coordinates on the ice sheet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indexICE = inpolygon(GRACEX,GRACEY,ice.MERRA.X,ice.MERRA.Y);
indexICEmerra = inpolygon(pM.X,pM.Y,ice.MERRA.X,ice.MERRA.Y);

indexICEnan = double(indexICE);
indexICEnan(indexICEnan==0) = nan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Index coordinates by subregion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subREG.index  = cell(1,4);
for i = 1:length(subREG.index)
    subREG.GRACE.index{i} = inpolygon(GRACEX,GRACEY,subREG.MERRA.XBUF{i},subREG.MERRA.YBUF{i});
    subREG.MERRA.index{i} = inpolygon(pM.X,pM.Y,subREG.MERRA.XBUF{i},subREG.MERRA.YBUF{i});
end

%% SAVE DATA
load(fullfile(datadir,'comparemassMerra'),'subREG','indexICE','indexICEmerra','indexICEnan',...
    'indexGL','indexGLnan','GRACEX','GRACEY','resamp','resamprate')

% for i = 1:4
%     subREG.MERRA.nanindex{i} = double(subREG.MERRA.index{i});
%     subREG.MERRA.nanindex{i}(~subREG.MERRA.nanindex{i}) = nan;
% end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE MELT DAYS
%   1) Near surface Air Temp >= 276.13 K
%   2)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freezing = 276.13;

T2M = load(fullfile(matDir,'T2M2003-2017'),'data','t','spacelim');
T2M.data = squeeze(T2M.data);
ST = load(fullfile(matDir,'TS2003-2017'),'data','t','spacelim');
ST.data = squeeze(ST.data);
QV2M = load(fullfile(matDir,'QV2M2003-2017'),'data','t','spacelim');
QV2M.data = squeeze(QV2M.data);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN PLOTTING:
%   1) comparison of entire GRACE time period
%   2)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
clf
subplot(1,2,1)
imPlot(sum(ST.data,3),[])
subplot(1,2,2)
imPlot(sum(meltData.allmeltMap,3),[])

%% SHOW MEAN CLIMATE VARIABLES VS MASS LOSS OVER GRACE
t2m = imgaussfilt(T2M.data,1);
st = imgaussfilt(ST.data,1);
qv2m = imgaussfilt(QV2M.data,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DIVIDE BY SEASON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
season = alldates;
season(month(alldates)>=12 | month(alldates)<=2) = 1;
season(month(alldates)>=3 & month(alldates)<=5) = 2;
season(month(alldates)>=6 & month(alldates)<=8) = 3;
season(month(alldates)>=9 & month(alldates)<=11) = 4;

t2m = t2m(:,:,season==3);
st = st(:,:,season==3);
qv2m = qv2m(:,:,season==3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get total maps for GRACE and MERRA (filtered)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% totalmelt = sum(t2mMelt>freezing,3);
meant2m = mean(t2m,3);
meanst = mean(st,3);
meanqv2m = mean(qv2m,3);

% totalmelt = imgaussfilt(sum(T2M.data>freezing,3),1);
% totalmelt = imgaussfilt(sum(meltData.allmeltMap,3),1);
GRACEdiff = massData.D(:,:,end) - ...
    massData.D(:,:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resample the maps for GRACE and MERRA (filtered)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mass = GRACEdiff(resamp,resamp);
mass = mass.*indexGL;

t2mmelt = (interp2(pM.X,pM.Y,meant2m',GRACEX,GRACEY));
t2mmelt = t2mmelt.*indexGL;

stmelt = (interp2(pM.X,pM.Y,meanst',GRACEX,GRACEY));
stmelt = stmelt.*indexGL;

qv2mmelt = (interp2(pM.X,pM.Y,meanqv2m',GRACEX,GRACEY));
qv2mmelt = qv2mmelt.*indexGL;

%% COMPARE TOTAL TEMP AND MASS FOR ALL OF GREENLAND

figure(2)
clf
suptitle('Mass loss to T2M comparison, 2003-2017')
region = {[1 4],[2 3]};
c = {'k','r'};
yl = [];

for j = 1:2
    if j == 1
        ax = subplot(2,1,1);
        hold on;grid on
        thisindex = indexICEnan;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TOTAL MASS LOSS v MEAN T2M
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        y = t2mmelt(:).*thisindex(:);
        x = -mass(:).*thisindex(:);
        %         h =  scatter(x,y,20,...
        %             distweight(:).*thisindex(:),'x','MarkerEdgeAlpha',0.2);
        % MAKE BINS
        binedge = -1500:250:5000;
        y2 = zeros(1,length(binedge)-1);
        x2 = zeros(1,length(binedge)-1);
        for i=1:length(binedge)-1
            y2(i) = nanmean(y(x>binedge(i) & x<binedge(i+1)));
            std2 = nanstd(y(x>binedge(i) & x<binedge(i+1)));
            x2(i) = binedge(i)+1/2*(binedge(i+1)-binedge(i));
            hebT2M = errorbar(x2(i),y2(i),std2,'o',...
                'color',[0.9 0 0],'linewidth',1,...
                'markerfacecolor',[0.9 0 0]);
        end
        xx2 = x(~isnan(y));
        yy2 = y(~isnan(y));
        [mT2M,f] = linear_m(xx2,yy2);
        hfitT2M = plot(xx2,f,'color',[0.9 0 0],'linewidth',1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TOTAL MASS LOSS v MEAN ST
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        y = stmelt(:).*thisindex(:);
        x = -mass(:).*thisindex(:);
        %         h =  scatter(x,y,20,...
        %             distweight(:).*thisindex(:),'x','MarkerEdgeAlpha',0.2);
        % MAKE BINS
        binedge = -1500:250:4000;
        y2 = zeros(1,length(binedge)-1);
        x2 = zeros(1,length(binedge)-1);
        for i=1:length(binedge)-1
            y2(i) = nanmean(y(x>binedge(i) & x<binedge(i+1)));
            std2 = nanstd(y(x>binedge(i) & x<binedge(i+1)));
            x2(i) = binedge(i)+1/2*(binedge(i+1)-binedge(i));
            hebST = errorbar(x2(i),y2(i),std2,'o',...
                'color',[0 0 0.6],'linewidth',1,...
                'markerfacecolor',[0 0 0.6]);
        end
        xx2 = x(~isnan(y));
        yy2 = y(~isnan(y));
        [mST,f] = linear_m(xx2,yy2);
        hfitST = plot(xx2,f,'color',[0 0 0.6],'linewidth',1);
        
        ylim([240 280]);
        xlim(minmax(binedge))

                ylabel('\circK')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TOTAL MASS LOSS v MEAN QV2M
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        yyaxis right
        y = qv2mmelt(:).*thisindex(:);
        x = -mass(:).*thisindex(:);
        %         h =  scatter(x,y,20,...
        %             distweight(:).*thisindex(:),'x','MarkerEdgeAlpha',0.2);
        % MAKE BINS
        binedge = -1500:250:5000;
        y2 = zeros(1,length(binedge)-1);
        x2 = zeros(1,length(binedge)-1);
        for i=1:length(binedge)-1
            y2(i) = nanmean(y(x>binedge(i) & x<binedge(i+1)));
            std2 = nanstd(y(x>binedge(i) & x<binedge(i+1)));
            x2(i) = binedge(i)+1/2*(binedge(i+1)-binedge(i));
            hebQV2M = errorbar(x2(i),y2(i),std2,'o',...
                'color',[0.9290 0.6940 0.1250],'linewidth',1,...
                'markerfacecolor',[0.9290 0.6940 0.1250]);
        end
        xx2 = x(~isnan(y));
        yy2 = y(~isnan(y));
        [mQV2M,f] = linear_m(xx2,yy2);
        hfitQV2M = plot(xx2,f,'color',[0.9290 0.6940 0.1250],'linewidth',1);
        ylabel('kg/kg')
        
        set(gca,'ylim',[-0.4 2.3]*1E-3)
                
        yl = [yl get(gca,'ylim')];
        
        colormap(ax,jet(100))
        xlabel('Kg per m^2 of mass loss')
        
        lgd = legend([hebT2M, hfitT2M,hebST, hfitST, hebQV2M, hfitQV2M],'mean T2M',...
            sprintf('Type 1 LSR (m = %0.2f\\circK/(1E3Kg/m^2); b = %0.2f)',mT2M(2)*1E3,mT2M(1)),...
            'mean ST',...
            sprintf('Type 1 LSR (m = %0.2f\\circK/(1E3Kg/m^2); b = %0.2f)',mST(2)*1E3,mST(1)),...
            'mean QV2M',...
            sprintf('Type 1 LSR (m = %0.e kg/kg/(1E3Kg/m^2); b = %0.2e)',mQV2M(2)*1E3,mQV2M(1)));

        lgd.Location = 'eastoutside';
                title(lgd,'Binned by mass loss (\DeltaM = 250 Kg/m^2; \bullet = \mu \pm\sigma)')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IMAGES OF DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax = subplot(2,2,j+2);
    hold on
    if j == 1
        C = mass.*indexICEnan;
        caxis([-4000 2000]);
        cmap = colormap(ax,bluewhitered(100,1));
    else
        C = stmelt.*indexICEnan;
    end
    Z = ones(size(C));
    surf(ptsGL.xp(resamp,resamp).*indexICEnan,ptsGL.yp(resamp,resamp).*indexICEnan,Z,C,...
        'edgecolor','none');
    cb = colorbar;
    
    if j==1;ylabel(cb,'Kg per m^2 of mass loss');end
    if j==2;ylabel(cb,'mean T2M\circK');end
    
    
    plotGLbackground( pG,ice,cmap,summit );
    axis image off
    set(gca,'ydir','reverse')
end
% ax = subplot(2,2,1:2);

%% COMPARE TOTAL TEMP AND MASS FOR DIFFERENT REGIONS OF GREENLAND

figure(2)
clf
suptitle('Mass loss to T2M comparison, 2003-2017')
region = {[1 4],[2 3]};
c = {'k','r'};
yl = [];

for j = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EAST AND WEST MASS LOSS v MELTDAY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax = subplot(2,2,j);
    hold on;grid on
    for reg = region{j}
        thisindex = subREG.GRACE.index{reg}.*indexICE;
        y = t2mmelt(:).*thisindex(:);
        x = -mass(:).*thisindex(:);
        %         h =  scatter(x,y,20,...
        %             distweight(:).*thisindex(:),'x','MarkerEdgeAlpha',0.2);
        % MAKE BINS
        %         binedge = 10.^(-2:0.5:3);
        binedge = -1000:250:5000;
        %         binedge = [200:300];
        yy2 = [];
        xx2 = [];
        for i=1:length(binedge)-1
            y2 = (y(x>binedge(i) & x<binedge(i+1)));
            yy2(i) = nanmean(y2);
            x2(i) = binedge(i)+1/2*(binedge(i+1)-binedge(i));
            errorbar(x2(i),yy2(i),nanstd(y2),'o',...
                'color',c{find(reg==region{j})},'linewidth',1,...
                'markerfacecolor',c{find(reg==region{j})})
        end
        if reg==2
            yy2 = yy2(~isnan(yy2));
            xx2 = x2(~isnan(yy2));
            yy2 = yy2(xx2<1000);
            xx2 = xx2(xx2<1000);
            
            [m,f] = linear_m(xx2,yy2);
            plot(xx2,f)
        end
        axis tight
        
        yl = [yl get(gca,'ylim')];
    end
    
    colormap(ax,jet(100))
    ylabel('mean T2M\circK')
    xlabel('Kg per m^2 of mass loss')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EAST AND WEST REGIONS AND DISTANCE COLORING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax = subplot(2,2,j+2);
    hold on
    if j == 1
        C = mass.*indexICEnan;
        caxis([-4000 2000]);
        cmap = colormap(ax,bluewhitered(100,1));
    else
        C = t2mmelt.*indexICEnan;
        %         caxis([0 3000]);
        %         cmap = colormap(ax,jet);
    end
    Z = ones(size(C));
    surf(ptsGL.xp(resamp,resamp).*indexICEnan,ptsGL.yp(resamp,resamp).*indexICEnan,Z,C,...
        'edgecolor','none');
    cb = colorbar;
    
    if j==1;ylabel(cb,'Kg per m^2 of mass loss');end
    if j==2;ylabel(cb,'mean T2M\circK');end
    
    
    plotGLbackground( pG,ice,cmap,summit );
    axis image off
    set(gca,'ydir','reverse')
    plotGLforeground(subREG, region{j},[],{'k','k','r','r'});
end
ax = subplot(2,2,1);
ylim([240 275]);
ax = subplot(2,2,2);
ylim([240 275]);

%%
yr = 2012;
thismeltdata = meltData.allmeltMap(:,:,monthnum(6,yr,meltData.alldates):...
    monthnum(10,yr,meltData.alldates));
meant2m = imgaussfilt(sum(thismeltdata,3),1);
GRACEdiff = massData.D(:,:,monthnum(10,yr,massData.thedates)) - ...
    massData.D(:,:,monthnum(6,yr,massData.thedates));
suptitle(sprintf('Mass loss to melt day comparison, JUN-SEPT %s',num2str(yr)))

%% ANOMALY
[ massAnom, massAVG ] = anomMonth(massData.D, massData.thedates,'GRACE');
figure(7)

clf
% cax = [0 1];

for m = 1:12
    ax = subplot(4,3,m);
    imj = massAVG(:,:,m);
    imagesc(imj(resamp,resamp).*indexICE)
    colormap(ax,bluewhitered(100,1))
    hold on
    %     plot(bx,by,'k--')
    %     plot(gx,gy,'k')
    title(datestr(datenum(1,m,1),'mmmm'))
    colorbar('eastoutside')
    axis tight
    axis off
end

y = unique(year(alldates));

suptitle(sprintf('Average Mass Change, kg per m^2 %i-%i',...
    min(y),max(y)))



%% SUMMER 2012

% mnth = [5 5];
% yr = [2003 2017];

mnth = [6 10];
yr = [2012 2012];

[meltJJA2012, massJJA2012] = meltdayvsmass(mnth,yr,massData.D(resamp,resamp,:),...
    massData.thedates,GRACEX,GRACEY,meltData.allmeltMap,meltData.alldates,pM.X,pM.Y);

[meltJJA2012anom, ~] = meltdayvsmass(mnth,yr,massAnom(resamp,resamp,:),...
    massData.thedates,GRACEX,GRACEY,meltData.normMeltMap,meltData.alldates,pM.X,pM.Y);


figure(5)
clf
suptitle(sprintf('%d/%d-%d/%d',mnth(1),yr(1),mnth(2),yr(2)))
ax3 = subplot(2,4,1);
hold on
fill(pG.gx./resamprate+0.5,pG.gy./resamprate+0.5,[0.6 0.6 0.6],'EdgeColor','none')
set(gca,'ydir','reverse')
imagesc(massJJA2012.*indexICE,'AlphaData',(indexICE));
colormap(ax3, bluewhitered([],1))
plot(pG.bx./resamprate+0.5,pG.by./resamprate+0.5,'--k','linewidth',0.2)
axis image off
title('Mass')

ax2 = subplot(2,4,2);
hold on
fill(pG.gx./resamprate+0.5,pG.gy./resamprate+0.5,[0.6 0.6 0.6],'EdgeColor','none')
set(gca,'ydir','reverse')
imagesc(meltJJA2012anom.*indexICE,'AlphaData',(indexICE));
colormap(ax2,parula)
plot(pG.bx./resamprate+0.5,pG.by./resamprate+0.5,'--k','linewidth',0.2)
axis image off
title('Melt Day Anomaly')

ax2 = subplot(2,4,5:6);
hold on
plotGLbackground( pG,ice,[0.7,1,1] );
plotGLforeground(subREG,1:4,1);
set(gca,'ydir','reverse')
axis image off
subplot(1,2,2)

[~,x,y] = cdfMelt(meltJJA2012, massJJA2012, indexICE, subREG );


% yyaxis right
%
% [a,i] = sort(meltJJA2012(:).*indexICE(:));
% b = distweight(:).*indexICEnan(:);
% b = b(i);
%
% plot(a(~isnan(b)),smooth(b(~isnan(b)),20))
% axis tight square
legend('nw','ne','se','sw','distance from edge (filt)','location','southoutside')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD STUFF
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISTANCE BETWEEN CENTER AND EDGE OF BUFFER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE A RELATIVE DISTANCE FROM EDGE OF ICESHEET TO SUMMIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance to summit
d2c = distance(resampleLAT,resampleLON,summitLONLAT(2),summitLONLAT(1));
% distance to edge of buffer
sz = size(GRACEX);
[~, d2b] = dsearchn([ice.MERRA.X,ice.MERRA.Y],[GRACEX(:),GRACEY(:)]);
d2b = reshape(d2b,sz);
% relative distance measure
distweight = d2b./(d2c + d2b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the relative distance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
clf
hold on
set(gca,'ydir','reverse')
cmap = colormap(jet);
plotGLbackground( pG,ice,cmap,summit );
surf(ptsGL.xp(resamp,resamp).*indexICEnan,ptsGL.yp(resamp,resamp).*indexICEnan,distweight.*indexICEnan);
colorbar
axis image off

%% SHOW TOTAL AIR TEMP MELT DAYS OVER GRACE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get total maps for GRACE and MERRA (filtered)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meant2m = imgaussfilt(sum(meltData.allmeltMap,3),1);
GRACEdiff = massData.D(:,:,end) - ...
    massData.D(:,:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resample the maps for GRACE and MERRA (filtered)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mass = GRACEdiff(resamp,resamp);
mass = mass.*indexGL;

t2mmelt = (interp2(pM.X,pM.Y,meant2m',GRACEX,GRACEY));
t2mmelt = t2mmelt.*indexGL;


%% COMPARE TOTAL TEMP AND MASS FOR DIFFERENT SIDES OF GREENLAND

figure(3)
clf
suptitle('Mass loss to melt day comparison, 2003-2017')
region = {[1 4],[2 3]};
c = {'k','r'};
yl = [];

for j = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EAST AND WEST MASS LOSS v MELTDAY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax = subplot(2,2,j);
    hold on;grid on
    for reg = region{j}
        thisindex = subREG.GRACE.index{reg}.*indexICE;
        x = t2mmelt(:).*thisindex(:);
        y = -mass(:).*thisindex(:);
        h =  scatter(x,y,20,...
            distweight(:).*thisindex(:),'x','MarkerEdgeAlpha',0.2)
        % MAKE BINS
        %         binedge = 10.^(-2:0.5:3);
        %         binedge = [0.1 1 5 10 50 100];
        binedge = [200:300];
        
        for i=1:length(binedge)-1
            y2 = (y(x>binedge(i) & x<binedge(i+1)));
            x2 = binedge(i)+1/2*(binedge(i+1)-binedge(i));
            %             errorbar(x2,nanmedian(y2),nanstd(y2),'o',...
            %                 'color',c{find(reg==region{j})},'linewidth',1,...
            %                 'markerfacecolor',c{find(reg==region{j})})
        end
        axis tight
        
        yl = [yl get(gca,'ylim')];
    end
    %     set(gca,'xscale','log')%,'xtick',10.^(-2:2),'xticklabels',...
    %         xticks)
    colormap(ax,jet(100))
    xlabel('# days T2M>0\circC, gaussian filter 1\sigma')
    ylabel('Kg per m^2 of mass loss')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EAST AND WEST REGIONS AND DISTANCE COLORING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,j+2)
    hold on
    Z = distweight;
    Z = Z.*indexICEnan;
    %     surf(ptsGL.xp(resamp,resamp).*indexICEnan,ptsGL.yp(resamp,resamp).*indexICEnan,Z);
    cmap = colormap(jet);
    plotGLbackground( pG,ice,cmap,summit );
    axis image off
    set(gca,'ydir','reverse')
    plotGLforeground(subREG, region{j},[],{'k','k','r','r'});
end
ax = subplot(2,2,1);
ylim(minmax(yl));
ax = subplot(2,2,2);
ylim(minmax(yl));
