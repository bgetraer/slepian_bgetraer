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

% Total cloud area fraction
CLDTOT = load(fullfile(matDir,'CLDTOT2003-2017'),'data','t','spacelim');



%% CROSS PLOTS 
plot(TsurfA)

%%
TA = meanTA(1:length(alldates));

[spa,~] = downsample_ts(SPA(:)./std(SPA),alldates);
[qva,~] = downsample_ts(QVA(:)./std(QVA),alldates);
[ta,t_ds] = downsample_ts(TA(:)./std(TA),alldates);

% Get 'thedates','ESTtotal','ESTtotalresid','total','alphavarall' from GREENLAND60.m
Gdata = load(fullfile(datadir,'Greenland60data'));
%%

massresid = interp1(Gdata.thedates,Gdata.ESTtotalresid,t_ds);
mass = interp1(Gdata.thedates,Gdata.total,t_ds);
mass = (mass - mean(mass));
[~,fl] = linear_m(t_ds,mass,2);
[massa,massm] = anomMonth(mass,t_ds);
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
season = t_ds;
season(month(t_ds)>=1 & month(t_ds)<=2) = 1;
season(month(t_ds)>=3 & month(t_ds)<=5) = 2;
season(month(t_ds)>=6 & month(t_ds)<=9) = 2;
season(month(t_ds)>=10 & month(t_ds)<=12) = 1;

theseason = 3;

figure(4)
clf
hold on
% plot(alldates,SPA./std(SPA))
% plot(alldates,QVA./std(QVA))
% plot(alldates,TA./std(TA))
% plot3(SPA./std(SPA),...
%     QVA./std(QVA),...
%     TA./std(TA),'.')


plot(t_ds(season~=theseason),ta(season~=theseason),'x')
plot(t_ds(season~=theseason),ma(season~=theseason),'x')


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