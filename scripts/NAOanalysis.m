NAO = dlmread(fullfile(datadir,'NAO_monthly.txt'));
NAOdates = datenum(NAO(:,1),NAO(:,2),repmat(15,size(NAO(:,1))));
NAOdata = NAO(:,3);

% do a 3 year filter, as pattern in residuals seems to be 3 year
windowSize = 3*12; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
NAO3yearfilt = filter(b,a,NAOdata);

% do a 6 month filter, as pattern in residuals seems to be 6 months
windowSize = 6; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
NAO6monthfilt = filter(b,a,NAOdata);
% do a centered moving average if you want it really really smooth
NAO6monthconv = conv(NAO6monthfilt, ones(1,3), 'same');

% crop by year
NAOt = NAOdates(year(NAOdates)>2002);
NAOy = NAO6monthfilt(year(NAOdates)>2002);
NAOy2 = NAO6monthfilt2(year(NAOdates)>2002);
NAOyconv = NAO6monthconv(year(NAOdates)>2002);

% rescale
NAOyprime = (NAOy-mean(NAOy))./std(NAOy);
NAOy2prime = (NAOy2-mean(NAOy2))./std(NAOy2);
NAOyconvprime = (NAOyconv-mean(NAOyconv))./std(NAOyconv);

%compare
GRACEresid = total(:)-f1all(:);
GRACEresidprime = (GRACEresid-mean(GRACEresid))/std(GRACEresid);


figure(1)
clf
hold on

ylimits = [-3 3];
ylim(ylimits);

yyaxis left
plot(NAOt,NAOyprime,'--','linewidth',1)
% plot(NAOt,NAOyconvprime);
ytickNAO = -0.6:0.2:0.4;
yticks = (ytickNAO-mean(NAOy))/std(NAOy);
set(gca,'ytick',yticks*3,'ylim',ylimits,'yticklabel',num2str(ytickNAO'),...
    'fontsize',12)
ylabel('North Atlantic Oscillation Index (m)','interpreter','latex')

yyaxis right
plot(thedates,GRACEresidprime,'linewidth',2)
ytickunit = [ylimits(1):ylimits(end)]./3;
ytickGRACE = round((ytickunit.*std(GRACEresid))+mean(GRACEresid));
set(gca,'ylim',ylimits,'yticklabel',num2str(ytickGRACE'),'fontsize',12)
ylabel('Residuals of 1$^{st}$ order polynomial regression (Gt)','interpreter','latex')

datetick
xlabel('Year','interpreter','latex')
title('\textbf{Correlation of NAO and Greenland mass signal}','interpreter','latex')
legend('NAO','GRACE residuals','location','northwest')