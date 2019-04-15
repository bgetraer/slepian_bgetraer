function [dataGRIS,meandataGRIS,dataGRISA,fp,fl,ml,x,y] = plotoverGRIS(time,data,mask )
%PLOTOVERGRIS Takes MERRA data over Greenland and generates basic plots of
%the data averaged over the extent of the icesheet.
%
% INPUT
%   time    the dates
%   data    the MERRA data matrix
%
% OUTPUT
%   dataGRIS    timeseries of the data averaged over the GRIS
%   meandataGRIS    the mean of the timeseries
%   dataGRISA       the seasonal anomaly (periodic signal removed)
%   fp          fitted seasonal signal 
%   fl          fitted linear signal
%   x,y         one cycle of the fitted seasonal signal
%
% Last Modified by bgetraer@princeton.edu 3/31/2019


% dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
% datadir = fullfile(dir,'datafiles/');
% load(fullfile(datadir,'comparemassMerra'),'indexICEmerra')
% mask = double(indexICEmerra);
% mask(~indexICEmerra)=nan;

% over the ice extent
% data = data.*mask';
%vectorize the data
datavec = reshape(data,size(data,1)*size(data,2),size(data,3));
% average over all of GRIS for each time step
dataGRIS = nanmean(datavec,1);
% mean and mean removed
meandataGRIS = nanmean(dataGRIS);
dataGRIS = dataGRIS-meandataGRIS;

% remove seasonal cycles
[~,fp] = periodic_m(time,dataGRIS,[365,365/2]);
dataGRISA = dataGRIS'-fp;

[ml,fl] = linear_m(time, dataGRISA);

x = 1:365;
y = fp(x);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   DIVIDE BY SEASON
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% season = time;
% season(month(time)>=1 & month(time)<=2) = 1;
% season(month(time)>=3 & month(time)<=5) = 2;
% season(month(time)>=6 & month(time)<=9) = 3;
% season(month(time)>=10 & month(time)<=12) = 4;
% theseason = 3;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   PLOT
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(3,1,1)
% hold on
% % plot(alldates(year(alldates)<2017), meanT(year(alldates)<2017)+meanTmean,'-')
% % scatter(alldates,SP+meanSP,1,season,'x')
% plot(time,dataGRIS+meandataGRIS)
% plot(time, fp+meandataGRIS,'-','linewidth',1)
% datetick
% grid on
% ylabel('kPa')
% axis tight
% 
% subplot(3,3,4:5)
% hold on
% cmap = repmat([0 0.4470 0.7410],4,1);
% cmap(theseason,:) = 0;
% colormap(cmap)
% scatter(time, dataGRISA,1,season,'x')
% % plot(time, dataGRISA,'-')
% % hf = plot(time, fl,'--r','linewidth',1);
% % legend(hf,sprintf('m = %0.2f\\circC per decade',ml(2).*365*10))
% datetick
% 
% grid on
% axis tight
% ylabel('kPa')
% 
% 
% subplot(3,3,6)
% hold on
% x = 1:365;
% y = fp(x);
% plot(x,y+meandataGRIS,'-','linewidth',1.5)
% grid on
% axis tight
% datetick('x','mmm')
% title('seasonal cycle removed')
% 
% 
% 
% 
% subplot(3,1,3)
% hold on
% yr = unique(year(time));
% BP = zeros(120,length(yr))+nan;
% 
% seasonmean = mean(dataGRISA(season==theseason)+meandataGRIS);
% for i = 1:length(yr)
%     yeardata = dataGRISA(season==theseason & year(time)==yr(i))+meandataGRIS;
%     BP(1:length(yeardata),i) = yeardata-seasonmean;
% end
% 
% boxplot(BP,yr,'symbol','rx')
% grid on
% title('JJAS summer anomalies')

end


