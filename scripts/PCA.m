%**************************************************************************
% Script for organizing all of the data I want to look at
%
% Last modified by: bgetraer@princeton.edu 4/2/2019
%**************************************************************************

% locate slepian_bgetraer function and datafile directories, and set workspace
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
datadir = fullfile(dir,'datafiles/');
merradir = fullfile(datadir,'MERRA2');
matdir = fullfile(merradir,'MerraMat');
addpath(dir,datadir);
addpath(fullfile(dir,'scripts/atmospheric_analysis/functions'))
addpath(fullfile(dir,'functions'))
setworkspace();

%% ORGANIZE MERRA DATA AVAILABLE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DAILY VALUES  2003 - 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface skin temperature 'TS' units K
TS = load(fullfile(matdir,'TS2003-2017'),'data','t','spacelim');

% Near surface 2m temperature 'T2M' units K
T2M = load(fullfile(matDir,'T2M2003-2017'),'data','t','spacelim');

% Temperature at 500 hPa pressure level 'T500' units K

% Specific humidity at 2m 'QV2M' units kg/kg
QV2M = load(fullfile(matDir,'QV2M2003-2017'),'data','t','spacelim');

% Specific humidity at 10m 'QV10M' units kg/kg    

% Specific humidity at 500 hPa 'Q500' units kg/kg

% Surface pressure 'PS' units kPa (I scaled the units)
PS = load(fullfile(matDir,'SP2003-2017'),'data','t','spacelim');

% 500 hPa pressure level height 'H500' units m

% Total cloud area fraction (1)
CLDTOT = load(fullfile(matDir,'CLDTOT2003-2017'),'data','t','spacelim');

% Surface net downward shortwave flux (W/m^2)
SWGNT = load(fullfile(matDir,'SWGNT2003-2017'),'data','t','spacelim');

% Maximum Precipitation Rate TPRECMAX

%% PREPARE MERRA DATA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Anomalies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
TSA_m = imgaussfilt(anomMonth(TS.data,alldates),2.5);
T2MA_m = imgaussfilt(anomMonth(T2M.data,alldates),2.5);
QV2MA_m = imgaussfilt(anomMonth(QV2M.data,alldates),2.5);
PSA_m = imgaussfilt(anomMonth(PS.data,alldates),2.5);
CLDTOTA_m = imgaussfilt(anomMonth(CLDTOT.data,alldates),2.5);
SWGNETA_m = imgaussfilt(anomMonth(SWGNT.data,alldates),2.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation onto GRACE grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TSA = zeros(128);
T2MA = zeros(128);
QV2MA = zeros(128);
PSA = zeros(128);
CLDTOTA = zeros(128);
SWGNETA = zeros(128);

for i = 1:length(alldates)
    TSA(:,:,i) = interp2(pM.X,pM.Y,TSA_m(:,:,i)',GRACEX,GRACEY);
    T2MA(:,:,i) = interp2(pM.X,pM.Y,T2MA_m(:,:,i)',GRACEX,GRACEY);
    QV2MA(:,:,i) = interp2(pM.X,pM.Y,QV2MA_m(:,:,i)',GRACEX,GRACEY);
    PSA(:,:,i) = interp2(pM.X,pM.Y,PSA_m(:,:,i)',GRACEX,GRACEY);
    CLDTOTA(:,:,i) = interp2(pM.X,pM.Y,CLDTOTA_m(:,:,i)',GRACEX,GRACEY);
    SWGNETA(:,:,i) = interp2(pM.X,pM.Y,SWGNETA_m(:,:,i)',GRACEX,GRACEY);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Downsample to month
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tsa,tds]  = downsample_ts(TSA, alldates);
t2ma        = downsample_ts(T2MA, alldates);
qv2ma       = downsample_ts(QV2MA, alldates);
psa         = downsample_ts(PSA, alldates);
cldtota     = downsample_ts(CLDTOTA, alldates);
swgneta     = downsample_ts(SWGNETA, alldates);

%% GET THE MASS

% GRACE mass data
massData = load('im_seqSH');
dmass = massData.D - massData.D(:,:,1);
massanom = anomMonth(massData.D,massData.thedates);

[ ~, mass_trend ] = linearMap( massanom, massData.thedates);

massa1 = massanom - mass_trend;

massa = interp1(Gdata.thedates,permute(massa1,[3,1,2]),tds);
massa = permute(massa,[2,3,1]);
massa = massa(resamp,resamp,:);
 
%% RESCALE EVERYTHING and MASK
mass =    massa./nanstd(massa(:)).*GRISnan;
ts =    tsa./nanstd(tsa(:)).*GRISnan;
t2m =    t2ma./nanstd(t2ma(:)).*GRISnan;
qv2m =   qv2ma./nanstd(qv2ma(:)).*GRISnan;
ps =    psa./nanstd(psa(:)).*GRISnan;
cldtot =    cldtota./nanstd(cldtota(:)).*GRISnan;
swgnet =   swgneta./nanstd(swgneta(:)).*GRISnan;



%% divide by season
season = tds;
season(month(tds)>=1 & month(tds)<=2) = 1;
season(month(tds)>=3 & month(tds)<=5) = 2;
season(month(tds)>=6 & month(tds)<=9) = 3;
season(month(tds)>=10 & month(tds)<=12) = 4;
theseason = 3;

mass_s =    mass(:,:,season==theseason);
ts_s =    ts(:,:,season==theseason);
t2m_s =    t2m(:,:,season==theseason);
qv2m_s =   qv2m(:,:,season==theseason);
ps_s =    ps(:,:,season==theseason);
cldtot_s =    cldtot(:,:,season==theseason);
swgnet_s =   swgnet(:,:,season==theseason);

%%


pcaMat = [  ...
    mass_s(:),...
    ts_s(:),...
    t2m_s(:), ...
    qv2m_s(:), ...
    ps_s(:), ...
    cldtot_s(:), ...
    swgnet_s(:)...
    ];

[COEFF, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(pcaMat);

%
figure(1) 
clf
hold on
INTERESTING = EXPLAINED(EXPLAINED>2);
hwidth = 0.4;
n = size(pcaMat,2);
colour = jet(n);
for i = 1:3
    yp = 0;
    yn=0;
    x1 = i-hwidth;
    x2 = i+hwidth;
    for j = 1:n
        if COEFF(j,i)>0
        hfill{j} = fill([x1,x1,x2,x2],[0,COEFF(j,i),COEFF(j,i),0]+yp,colour(j,:));
        text(i+0.42,COEFF(j,i)/2+yp,sprintf('%0.3f',COEFF(j,i)))
        yp = COEFF(j,i)+yp;
        elseif COEFF(j,i)<0
        hfill{j} = fill([x1,x1,x2,x2],[0,COEFF(j,i),COEFF(j,i),0]+yn,colour(j,:));
        text(i+0.42,COEFF(j,i)/2+yn,sprintf('%0.3f',COEFF(j,i)))
        yn = COEFF(j,i)+yn;
        end
    end
end
% ylim([-1,3])
% xlim([0.4,2.6])
refline(0,0)
xlab = 'PC%0.0i: %0.2f%%';
set(gca,'YGrid', 'on','xtick',1:3,'xticklabel',...
    {sprintf(xlab,1,INTERESTING(1)),sprintf(xlab,2,INTERESTING(2)),...
    sprintf(xlab,3,INTERESTING(3))})
lgd = legend([hfill{1},hfill{2},hfill{3},hfill{4}],...hfill{5},hfill{6},hfill{7}],...
    'MASS', 'ST','T2M','QV2M','PS','CLDTOT','SWGNET');
lgd.Location = 'southwest';
%%
% Get 'thedates','ESTtotal','ESTtotalresid','total','alphavarall' from GREENLAND60.m
% Gdata = load(fullfile(datadir,'Greenland60data'));
% massresid = interp1(Gdata.thedates,Gdata.ESTtotalresid,t_ds);
% mass = interp1(Gdata.thedates,Gdata.total,t_ds);
% mass = (mass - mean(mass));
% [~,fl] = linear_m(t_ds,mass,2);
% [massa,massm] = anomMonth(mass,t_ds);
figure(3)
clf
hold on
% plot(t_ds,massa)

% plot(t_ds,fl)
plot(t_ds,massresid)
% plot(t_ds,massmodel)
MA = massresid;
ma = MA(:)./std(MA);
% plot(t_ds,MA)
% plot(t_ds,ma)
datetick

%%

figure(4)
clf
hold on
% plot(alldates,SPA./std(SPA))
% plot(alldates,QVA./std(QVA))
% plot(alldates,TA./std(TA))
% plot3(SPA./std(SPA),...
%     QVA./std(QVA),...
%     TA./std(TA),'.')

% 
% scatter(tds(season==theseason),permute(tsa(:,:,season==theseason),[3,1,2]),'x')
% scatter(tds(season==theseason),permute(massa(:,:,season==theseason),[3,1,2]),'x')
% datetick

%%
figure(4)
clf
hold on

mas = ma(season~=theseason);
tas = ta(season~=theseason);
spas = spa(season~=theseason);
qvas = qva(season~=theseason);

theseason = 1;
mas = ma(season==theseason);
tas = ta(season==theseason);
spas = spa(season==theseason);
qvas = qva(season==theseason);

pcaMat = [SPA(:)./std(SPA),QVA(:)./std(QVA),TA(:)./std(TA)];
pcaMat = [ma(:),spa(:),qva(:),ta(:)];
pcaMat = [mas(:),spas(:),qvas(:),tas(:)];
[COEFF, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(pcaMat);

INTERESTING = EXPLAINED(EXPLAINED>20);
hwidth = 0.4;
colour = jet(4);
for i = 1:length(INTERESTING)
    yp = 0;
    yn=0;
    x1 = i-hwidth;
    x2 = i+hwidth;
    for j = 1:4
        if COEFF(j,i)>0
        fill([x1,x1,x2,x2],[0,COEFF(j,i),COEFF(j,i),0]+yp,colour(j,:))
        text(i+0.42,COEFF(j,i)/2+yp,sprintf('%0.3f',COEFF(j,i)))
        yp = COEFF(j,i)+yp;
        elseif COEFF(j,i)<0
        fill([x1,x1,x2,x2],[0,COEFF(j,i),COEFF(j,i),0]+yn,colour(j,:))
        text(i+0.42,COEFF(j,i)/2+yn,sprintf('%0.3f',COEFF(j,i)))
        yn = COEFF(j,i)+yn;
        end
    end
end
ylim([-1,2])
xlim([0.4,2.6])
refline(0,0)
xlab = 'PC%0.0i: %0.2f%%';
set(gca,'YGrid', 'on','xtick',1:2,'xticklabel',...
    {sprintf(xlab,1,INTERESTING(1)),sprintf(xlab,2,INTERESTING(2))})
lgd = legend('MASS ANOMALY','SURFACE PRESSURE ANOMALY','500 hPa WATER VAPOR ANOMALY',...
    'ST ANOMALY');
lgd.Location = 'southwest';