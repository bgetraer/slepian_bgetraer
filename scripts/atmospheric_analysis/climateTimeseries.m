% locate slepian_bgetraer function and datafile directories, and set workspace
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
merra_dir = fullfile(dir,'datafiles/MERRA2');
addpath(fullfile(dir,'scripts/atmospheric_analysis/functions'))
addpath(fullfile(dir,'functions'))
setworkspace('/Users/benjamingetraer/Documents/IndependentWork/SH_Workspace');
datadir = fullfile(dir,'datafiles/');
matdir = fullfile(merra_dir,'MerraMat');


yrs = repmat((1980:2017),12,1);
months = repmat((1:12)',size(unique(yrs),1),1);
day = repmat(15,size(months,1),1);

allMdates = datenum(yrs(:),months,day);

% startdate = datenum(2003,01,01);
% enddate = datenum(2003,1,31);
% alldates = startdate:enddate;

%% LOAD MONTHLY DATA
% ncdfdir = fullfile(merra_dir,'GrnlandMnthly');
% [ thevardata, thetimedatastr, thespacelim, thelevdata ] = ...
%     getMerra2VAR(alldates,'T2MMIN',ncdfdir,'statM_2d_slv_Nx');

Temp = load(fullfile(matdir,'monthlydata'),'T2MMIN','T2MMAX','T2MMEAN');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Mean Temp of the Northern Hemisphere and Globe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NHmeancsv = csvread(fullfile(datadir,'NHemisMean.csv'),3);
NHmean = NHmeancsv(:,2:13);
NHmths = repmat(1:12,size(NHmean,1),1);
NHyr = repmat(NHmeancsv(:,1),1,size(NHmths,2));
NHdates = datenum(NHyr,NHmths,15);
[NHd,sortI] = sort(NHdates(:));
NHm = NHmean(sortI);
[NHm,NHd] = downsample_ts(NHm,NHd);

Glmeancsv = csvread(fullfile(datadir,'GlobalMean.csv'),3);
Glmean = Glmeancsv(:,2:13);
Glmths = repmat(1:12,size(Glmean,1),1);
Glyr = repmat(Glmeancsv(:,1),1,size(Glmths,2));
Gldates = datenum(Glyr,Glmths,15);
[Gld,sortI] = sort(Gldates(:));
Glm = Glmean(sortI);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the ENSO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ENSOcsv = textscan(fopen(fullfile(datadir,'ENSO.txt')),repmat('%s',1,13),'HeaderLines',4);
for i = 1:length(ENSOcsv)
    ENSOmat(:,i) = str2double(string(ENSOcsv{i}));
end
ENSO = ENSOmat(:,2:13);
ENSOmths = repmat(1:12,size(ENSO,1),1);
ENSOyr = repmat(ENSOmat(:,1),1,size(ENSOmths,2));
ENSOdates = datenum(ENSOyr,ENSOmths,15);
[ENSOd,ENSOI] = sort(ENSOdates(:));
ENSOsort = ENSO(ENSOI);

%% PLOT NH V GRIS ALL SEASONS

figure(1)
clf
% subplot(1,5,1:4)
hold on

% PLOT THE MEAN T2M OVER GREENLAND
tdata = squeeze(Temp.T2MMEAN).*GRISmerranan';
tdatam = squeeze(nanmean(nanmean(tdata)));
[tdatama, tdatamonth] = anomMonth(tdatam,allMdates);

hT6 = plot(allMdates,movmean(tdatama,6),'color',[0.8,0,0]);
t10yM = movmean(tdatama,120);
t10ySTD = movstd(tdatama,120);
[y,x]=downsample_ts(t10yM,allMdates,'monthly');

hTf = fill([allMdates; flip(allMdates)],[t10yM+t10ySTD; flip(t10yM-t10ySTD)],...
    'r','facealpha',0.1,'EdgeColor','none');
hT10 = plot(x,y,'color',[0.8,0,0],'linewidth',2);
datetick
grid on
set(gca,'yminorgrid','on')


% PLOT THE NH TEMP

ind = year(NHd) >= 1980 & year(NHd) < 2018;
NHd1980 = NHd(ind);
NHmstd = std(NHm(ind));
NHm1980 = (NHm(ind)-mean(NHm(ind)));

NH10yM = movmean(NHm1980,12*10);
NH10ySTD = movstd(NHm1980,12*10);

hNHf = fill([NHd1980; flip(NHd1980)],[NH10yM+NH10ySTD flip(NH10yM-NH10ySTD)],...
    'b','facealpha',0.1,'EdgeColor','none');

hNH10 = plot(NHd1980,NH10yM,'-','color',[0,0,0.8],'linewidth',2);
hNH6 = plot(NHd1980,movmean(NHm1980,6),'-','color',[0,0,0.8]);


ind = year(NHd1980)>=2012 & year(NHd1980)<2017;
[mlNH, flNH] = linear_m(NHd1980(ind),NHm1980(ind));

m = linear_m(NHd1980,NH10yM);
[r] = corrcoef(NHd1980,NH10yM);
r(2)^2
% refline(m(2),m(1))




datetick
xlim([datenum(1980,1,1),datenum(2017,12,31)])
ylabel('\circK')
title('Near Surface Temperature Anomaly, 1980-2017 (all seasons)')

legend([hT6,hT10,hTf,hNH6,hNH10,hNHf],...
    'MERRA T2M (monthly mean, 6m moving window)',...
    'MERRA T2M (10y moving mean)',...
    'MERRA T2M (10y moving std)',...
    'GISTEMP LOTI (monthly mean, 6m moving window)',...
    sprintf('GISTEMP LOTI (10y moving mean, slope \\approx %0.2f\\circK/decade)',m(2)*365*10),...
    'GISTEMP LOTI (10y moving std)','location','southoutside')

% % ENSO
% ind = year(ENSOd) >= 1980 & year(ENSOd) < 2018;
% ENSOp = ENSOsort;
% ENSOp(ENSOsort<=0) = 0;
% ENSOn = ENSOsort;
% ENSOn(ENSOsort>0) = 0;
% area(ENSOd(ind),ENSOp(ind),'facecolor',[0.8 0 0],'facealpha',0.3)
% area(ENSOd(ind),ENSOn(ind),'facecolor',[0 0 0.8],'facealpha',0.3)

% ax = subplot(1,5,5);
% % ax.Position = [0.7813    0.1100    0.1237    0.8150];
% % ax.Position = [0.1300  0.3587 0.6122 0.5663];
% ax.Position = [0.7813    0.3587  0.1237 0.5663];
% d
% hold on
% xlim([-4 4])
% histogram(tdatama(ind80s),15)
% histogram(tdatama(ind10s),15)
% view(90,-90)
%% mean in 80s and 10s
ind80s = year(allMdates)>=1980 & year(allMdates)<1995;
ind10s = year(allMdates)>=2003 & year(allMdates)<=2017;

ra = rand(100,1);
rb = rand(100,1);

ra = ra+5;
rb = rb-10;

figure(4)
clf
% mean(tdatama(ind80s))
% std(tdatama(ind80s))
% mean(tdatama(ind10s))
% std(tdatama(ind10s))
% [h,p,ci,stats] = ttest(tdatama(ind80s)-mean(tdatama(ind80s)))
mean(tdatama(ind80s)) - mean(tdatama(ind10s))
[h,p,ci,stats] = ttest2(tdatama(ind80s),tdatama(ind10s))
hold on
hist(tdatama(ind80s))
hist(tdatama(ind10s))
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DIVIDE BY SEASON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
season = allMdates;
season(month(allMdates)>=12 | month(allMdates)<=2) = 1;
season(month(allMdates)>=3 & month(allMdates)<=5) = 2;
season(month(allMdates)>=6 & month(allMdates)<=8) = 3;
season(month(allMdates)>=9 & month(allMdates)<=11) = 4;
theseason = 3;
%% PLOT NH V GRIS SUMMER


figure(1)
clf
hold on

% PLOT THE MEAN T2M OVER GREENLAND
tdata = squeeze(Temp.T2MMEAN).*GRISnan';
tdatam = squeeze(nanmean(nanmean(tdata)));
[tdatama, tdatamonth] = anomMonth(tdatam,allMdates);

color = [0,0.8,0.8;...
    0,0.8,0;...
    0.9,0.5,0.3;...
    0.7,0,0.7];

seas = {'Winter (DJF)','Spring (MAM)','Summer (JJA)','Fall (SON)'};

for i = 1:4
    subplot(4,1,i)
    hold on
    
    % deal with Dec Jan Feb year overlaps
    if i == 1
        x1 = allMdates(season==i);
        y1 = tdatama(season==i);
        
        for j = 1:length(x1)
            if month(x1(j))==12
                x1(j) = datenum(year(x1(j))+1,month(x1(j)),15);
            end
        end
        
        [x1,si] = sort(x1);
        y1 = y1(si);
        
        [y,x]=downsample_ts(y1,x1,'yearly');
        [stdy]=downsample_ts(y1,x1,'yearly','function','std');
        
        for j = 1:length(x)
                x(j) = datenum(year(x(j)),1,15);
        end
        
        x = x(1:end-1);
        y = y(1:end-1);
        stdy = stdy(1:end-1);
    else
        [y,x]=downsample_ts(tdatama(season==i),allMdates(season==i),'yearly');
        [stdy]=downsample_ts(tdatama(season==i),allMdates(season==i),'yearly','function','std');
    end
    
    hT10 = plot(x,y,'-o','color',color(i,:),'markeredgecolor','none',...
        'markerfacecolor',color(i,:),'linewidth',2);
    
    xx = allMdates(season==i);
    yy = tdatama(season==i);
    
%     plot(xx,yy,'x')
    
    [yS, iS] = sort(y,'descend');

    % # of years in top ten since 2003
%     sum(year(x(iS(1:10)))>=2003)
    
    year(x(iS(6:10)))
    yS(6:10)'
    stdy(iS(6:10))


%     mu = mean(yy);
%     sigma = std(yy);
%     n = length(yy);
%     N = 10000;
%     mTests = zeros(1,N);
%     
%     for j=1:N
%         Y = randn(n,1).*sigma + mu;
%         [mTest] = linear_m(xx,Y);
%         mTests(j) = mTest(2);
%     end
%     
    [mT,fT] = linear_m(xx,yy);
%     
%     yyy = yy - fT;
%     
%     mu = mean(yyy);
%     sigma = std(yyy);
%     n = length(yyy);
%     N = 10000;
%     ciTests = zeros(1,N);
%     
%     for j=1:N
%         Y = randn(n,1).*sigma + mu;
%         [ciTest] = linear_m(xx,Y);
%         ciTests(j) = ciTest(2);
%     end
%     % Compute centile
%     nless = sum(abs(mTests) < abs(mT(2)));
%     nequal = sum(abs(mTests) == abs(mT(2)));
%     centile(i) = 100 * (nless + 0.5*nequal) / length(mTests);
%     % compute confidence interval
%     coni95(i) = prctile(abs(ciTests),95);
%     coni90(i) = prctile(abs(ciTests),90);
    
    rT = corrcoef(x,y);
    totslope = refline(mT(2),mT(1));
    set(totslope,'color','k','linestyle','-','linewidth',1);
    
    [m,f] = linear_m(xx(year(xx)>=2003),yy(year(xx)>2002));
    rg = corrcoef(xx(year(xx)>=2003),yy(year(xx)>2002));
    graceslope = plot(xx(year(xx)>=2003),f,'--k','linewidth',1);
    
   
    xlim([datenum(1980,1,1),datenum(2017,12,31)])
    datetick('keeplimits')
    grid on
    set(gca,'yminorgrid','on')
    ylabel('\circK')
    if i==1
        legend([totslope,graceslope],...
            sprintf('LSR 1980-2017: m = %0.2f\\circK/decade R^2 = %0.2f',...
            mT(2)*365*10,rT(2)^2),...
            sprintf('LSR 2003-2017: m = %0.2f\\circK/decade R^2 = %0.2f',...
            m(2)*365*10,rg(2)^2),...
            'location','northwest')
    else
        legend([totslope,graceslope],...
            sprintf('m = %0.2f\\circK/decade R^2 = %0.2f',...
            mT(2)*365*10,rT(2)^2),...
            sprintf('m = %0.2f\\circK/decade R^2 = %0.2f',...
            m(2)*365*10,rg(2)^2),...
            'location','northwest')
    end
    title(seas{i})
end


%% GRACE
% % Get 'thedates','ESTtotal','ESTtotalresid','total','alphavarall' from GREENLAND60.m
% Gdata = load(fullfile(datadir,'Greenland60data'));
% GraceX = Gdata.thedates;
% GraceY = Gdata.ESTtotalresid;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   DIVIDE BY SEASON
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% theseason = 3;
%
% figure(5)
% clf
% hold on
% % plot(GraceX,GraceY,'linewidth',2);
% for i = 3
%     season = GraceX;
% season(month(GraceX)>=12 | month(GraceX)<=2) = 1;
% season(month(GraceX)>=3 & month(GraceX)<=5) = 2;
% season(month(GraceX)>=6 & month(GraceX)<=8) = 3;
% season(month(GraceX)>=9 & month(GraceX)<=11) = 4;
%     [y,x]=downsample_ts(GraceY(season==i),GraceX(season==i),'yearly');
%     hold on
%     yyaxis left
%     plot(x,y,'-','color',color(i,:),'linewidth',2);
%
%     season = alldates;
% season(month(alldates)>=12 | month(alldates)<=2) = 1;
% season(month(alldates)>=3 & month(alldates)<=5) = 2;
% season(month(alldates)>=6 & month(alldates)<=8) = 3;
% season(month(alldates)>=9 & month(alldates)<=11) = 4;
%
% yyaxis right
%     [y,x]=downsample_ts(tdatama(season==i),alldates(season==i),'yearly');
%     hT10 = plot(x(year(x)>2003),y(year(x)>2003),'--o','color',color(i,:),'markeredgecolor','none',...
%         'markerfacecolor',color(i,:),'linewidth',2);
%     grid on
%
%     % lgd = legend('ENSO+','ENSO-','N. Hemis. mean ST','GRIS mean ST',...
%     %     sprintf('GRIS ST 2003-17: %0.2f\\circC per decade',ml(2).*365.*10.*std(TSA)),...
%     % sprintf('N. Hemis. ST 2012-16: %0.2f\\circC per decade',mlNH(2).*365.*10.*NHmstd));
%     % lgd.Location = 'northwest';
%
%     ylabel('Gt')
%     axis tight
%     datetick
%     title('MEAN MONTHLY ANOMALY')
% end

%% COMPARE TEMP AND NAO

NAO_daily = load(fullfile(datadir,'NAO','NAO_daily.ascii'));
NAOdates = datenum(NAO_daily(:,1:3));
NAOdata = NAO_daily(:,4);

yyyy = 2012;
in2012 = year(NAOdates)>=2003 & year(NAOdates)<=2017; %& month(NAOdates)>=5 & month(NAOdates)<=9;
NAOx = NAOdates(in2012);
NAOy = NAOdata(in2012);


startdate = datenum(2003,01,01);
enddate = datenum(2017,12,31);
alldates = startdate:enddate;

[~,~,TSA] = ...
    plotoverGRIS(alldates,T2M.data.*GRISnan');
in2012 = year(alldates)>=2003 & year(alldates)<=2017; %&  month(alldates)>=5 & month(alldates)<=9;
Ty = TSA(in2012);

figure(7)
clf

subplot(1,3,1)
hold on
% plot(NAOy,Ty,'x')
naoy = (NAOy-nanmean(NAOy))/nanstd(NAOy);
ty = (Ty-nanmean(Ty))/nanstd(Ty);
histogram2(naoy,ty,50,'DisplayStyle','tile','FaceColor','flat')
[mlin, flin,DELTA] = linear_m(naoy(~isnan(naoy)),ty((~isnan(naoy))));
% plot(naoy(~isnan(naoy)),flin,'r','linewidth',2)
mline = refline(mlin(2),mlin(1));
set(mline,'color','r','linewidth',2)
r = corrcoef(naoy(~isnan(naoy)),ty((~isnan(naoy))));
rs = r(2)^2;

[a, ~, ~, ~, e]=pca([naoy(:) ty(:)]);
mpca = a(1,1)/a(2,1);
pcline = refline(mpca,0);
set(pcline,'color','b','linewidth',2)

legend([mline,pcline],...
    sprintf('m = %0.2f; R^2 = %0.2f',mlin(2),rs),...
    sprintf('m = %0.2f; explained var = %0.2f',mpca,e(1)))

xlabel('Daily NAO, 2003-2017, normalized')
ylabel('Daily T2M anomaly, 2003-2017, normalized')
% datetick
%%

P = load(fullfile(matDir,'SP2003-2017'),'data','t','spacelim');
P.data = squeeze(P.data)./1000;

[~,~,SPA] = ...
    plotoverGRIS(alldates,P.data.*GRISnan');
%% NAO vs a specific year

theyrs = [2012 2015];

figure(7)

for i = 1:2
    subplot(1,3,1+i)
    theyr = theyrs(i);
    yyaxis left
    cla
    hold on
%     
%     in2012 = year(alldates)>=theyr & year(alldates)<=theyr &  month(alldates)>=5 & month(alldates)<=9;
%     Py = SPA(in2012)./std(SPA(in2012));
    
    in2012 = year(NAOdates)>=theyr & year(NAOdates)<=theyr & month(NAOdates)>=5 & month(NAOdates)<=9;
    NAOx = NAOdates(in2012);
    NAOy = NAOdata(in2012);
    
    plot(NAOx,movmean(NAOy,1),'linewidth',1)
    ylabel('Daily NAO index')
    
    in2012 = year(alldates)>=theyr & year(alldates)<=theyr &  month(alldates)>=5 & month(alldates)<=9;
    Ty = TSA(in2012)./std(TSA(in2012));
    
    yyaxis right
    cla
    plot(alldates(in2012),movmean(Ty,1),'linewidth',1)
    set(gca,'ydir','reverse')
    ylabel('Daily 2M Air Temperature Anomaly \circK')
    datetick
    title(theyr)
end