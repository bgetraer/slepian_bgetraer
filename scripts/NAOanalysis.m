

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

clear year
% crop by year
NAOt = NAOdates(year(NAOdates)>2002);
NAOy = NAO6monthfilt(year(NAOdates)>2002);
NAOy3 = NAO3yearfilt(year(NAOdates)>2002);
NAOy2 = NAO6monthfilt2(year(NAOdates)>2002);
NAOyconv = NAO6monthconv(year(NAOdates)>2002);

% rescale
NAOyprime = (NAOy-mean(NAOy))./std(NAOy);
NAOy2prime = (NAOy2-mean(NAOy2))./std(NAOy2);
NAOy3prime = (NAOy3-mean(NAOy3))./std(NAOy3);
NAOyconvprime = (NAOyconv-mean(NAOyconv))./std(NAOyconv);

%compare
GRACEm1resid = total(:)-f1all(:);
GRACEm1residprime = (GRACEm1resid-mean(GRACEm1resid))/std(GRACEm1resid);
GRACEfresid = total(:)-ESTtotal(:);
GRACEfresidprime = (GRACEfresid-mean(GRACEfresid))/std(GRACEfresid);

%%

mainp = [0.1300    0.1100    0.7750    0.8150];
lowp = [0.1300    0.1100    0.7750    0.3];

f=figure(1);
clf

ylimitsmain = [-5.6 3];
yticklabelsNAOmain = round(-0.5:0.125:0.5,2);
yticklabelsGRACEmain = [-150:25:150];
ytickNAOmain = (yticklabelsNAOmain-mean(NAOy))/std(NAOy);
ytickGRACEmain = (yticklabelsGRACEmain-mean(GRACEm1resid))/std(GRACEm1resid);
%LEFT%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis left

nao1h = plot(NAOt,NAOyprime,'--','linewidth',1);
ylim(ylimitsmain);
xlim(minmax(NAOt));

set(gca,'ytick',ytickNAOmain*3,'ylim',ylimitsmain,'yticklabel',num2str(yticklabelsNAOmain'),...
    'fontsize',12,'tickdir','out')

ylabel('North Atlantic Oscillation Index (m)','interpreter','latex')

%RIGHT%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis right 

m1h = plot(thedates,GRACEm1residprime,'linewidth',2);

ylim(ylimitsmain);

set(gca,'ytick',ytickGRACEmain*3,'ylim',ylimitsmain,'yticklabel',num2str(yticklabelsGRACEmain'),...
    'fontsize',12)
ylabel('GRACE residuals (Gt)','interpreter','latex')
datetick
xlim(minmax(NAOt));

xlabel('Year','interpreter','latex')
grid on
title('\textbf{Patterns of NAO and Greenland mass signals}','interpreter','latex')


% lower axis
lowaxp = [0.1300    0.1100    0.7750    0.3];
lowax = axes('Parent',f,'Position',lowaxp);
nao3h = plot(NAOt,NAOy3prime,'-k','linewidth',1);
datetick
ylim([-3 3]);
ytickNAO3 = round(-0.3:0.1:0.1,2);
yticks = (ytickNAO3-mean(NAOy3))/std(NAOy3);
set(lowax,'xticklabels','','Color','none','box','off','yaxislocation',...
    'left','YColor','k','fontsize',12,'ytick',yticks*3,'yticklabels',...
    num2str(ytickNAO3'),'tickdir','out')

ylabel('','interpreter','latex')
xlim(minmax(NAOt));
grid on

lgd = legend([m1h nao1h nao3h],'$\underline{m}_1$','NAO 6-month moving window','NAO 3-year moving window','location','northwest');
lgd.Interpreter = 'latex';