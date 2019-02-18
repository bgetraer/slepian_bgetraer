NAO_daily = load('NAO_daily.ascii');
NAOdates = datenum(NAO_daily(:,1:3));
NAOdata = NAO_daily(:,4);

in2012 = year(NAOdates)==2012 & month(NAOdates)>=5 & month(NAOdates)<=9;

figure(10)
clf
% plot(NAOdates(in2012),NAOdata(in2012));
hold on


x = NAOdates(in2012);
y = NAOdata(in2012);
yM = movmean(y,5);

plot(x,yM)

[pkP,lcP] = findpeaks(yM,x,'minpeakprominence',0.25);
[pkN,lcN] = findpeaks(-yM,x,'minpeakprominence',0.25);

pk = [pkP' -pkN'];
lc = [lcP' lcN'];



plot(lcP, pkP+0.03,'rv','markerfacecolor','r')
plot(lcN, -pkN-0.03,'b^','markerfacecolor','b')

plot(datenum(2012,7,8),0,'o')

datetick

yM = movmean(NAOdata,5);
x = NAOdates;

%%
[pkP,lcP] = findpeaks(yM,x,'minpeakprominence',0.4);
[pkN,lcN] = findpeaks(-yM,x,'minpeakprominence',0.4);

NAOdata = [pkP' -pkN'];
NAOdates = [lcP' lcN'];

stdC = 2;
[x_P,y_P] = selectNAO(NAOdates,NAOdata,'summer','outlier1.5P','');
[x_N,y_N] = selectNAO(NAOdates,NAOdata,'summer','outlier1.5N','');
[x_A,y_A] = selectNAO(NAOdates,NAOdata,'summer',sprintf('outlier%fA',stdC),'');
[x_Ag,y_Ag] = selectNAO(NAOdates,NAOdata,'summer',sprintf('outlier%fA',stdC),'grace');

[x_S,y_S] = selectNAO(NAOdates,NAOdata,'summer','','');

edges = datenum(strcat('1/1/',num2str([1940:10:2030]')));

figure(2)
clf

yyaxis right
hold on
hbarP = histBAR(histcounts(x_P,edges),edges,'b');
hbarN = histBAR(histcounts(x_N,edges),edges,'red',-1);

yyaxis left
hold on

pN = plot(x_N,y_N,'ko','markerfacecolor','k');
pP = plot(x_P,y_P,'ko','markerfacecolor','k');
pA = plot(x_A,y_A,'ko','markersize',12);
pAg = plot(x_Ag,y_Ag,'rv','markersize',12);
[sigma,mu,pS] = stdline(x_S,y_S);

datetick

% legend
lgd = legend([pP,pA,pAg,pS,hbarP,hbarN],...
    {'Summer NAO outliers (\sigma)',...
    sprintf('Summer NAO outliers (%s\\sigma)',num2str(stdC)),...
    sprintf('Summer NAO outliers during GRACE (%s\\sigma)',num2str(stdC)),...
    sprintf('std(Summer NAO) (\\sigma = %0.2f, \\mu = %0.2f)',sigma, mu),...
    'Positive summer NAO outliers (\sigma)',...
    'Negative summer NAO outliers (\sigma)'});
lgd.Position = [0.5071 0.8214 0.3714 0.1726];

% center axes
yyaxis left
YL = get(gca, 'YLim');
maxlim = max(abs(YL));
set(gca, 'YLim', [-maxlim maxlim]);

yyaxis right
yR = gca;
maxlim = max(abs(yR.YLim));
yR.YLim = [-maxlim maxlim];
yR.YTickLabel = abs(yR.YTick);


%%
