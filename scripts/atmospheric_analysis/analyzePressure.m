%**************************************************************************
% Script for analyzing Pressure over the Greenland Ice Sheet
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

P = load(fullfile(matDir,'SP2003-2017'),'data','t','spacelim');
Pdata = P.data;
% P = save(fullfile(matDir,'SP2003-2017'),'data','t','spacelim');

%%
startdate = datenum(2003,01,01);
enddate = datenum(2017,12,31);
alldates = startdate:enddate;


figure(3)
clf
[SP,meanSP,SPA,fp,fl,ml,x,y] = plotoverGRIS(alldates,Pdata );


[ml1,fl1,~,~,t] = linear_m(alldates(year(alldates)>2008),...
    SPA(year(alldates)>2008));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,1)
hold on
ylabel('kPa')
title('mean surface P over GRIS')

subplot(3,3,4:5)
hold on
% plot(t, fl1,'-r','linewidth',2)
ylabel('kPa')
title('mean surface P anomaly over GRIS')


%% NORMALIZE PIXEL BY PIXEL

[ Tslopes ] = linearMap( data, alldates );

figure(5)
clf
imj = Tslopes;
imPlot(imj.*365)
colorbar
axis square off



%% MOVIE
% The GRIS Merra nan mask
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
datadir = fullfile(dir,'datafiles/');
load(fullfile(datadir,'comparemassMerra'),'indexICEmerra')
GRISnan = double(indexICEmerra);
GRISnan(~indexICEmerra)=nan;

figure(6)
clf
for i = 1:length(alldates)
imPlot(P.Pdata(:,:,i).*GRISnan')
caxis([62.9636  105.4815])
title(datestr(alldates(i)))
pause(0.01)
end