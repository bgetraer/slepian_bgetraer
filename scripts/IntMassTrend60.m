%**************************************************************************
% Integrated Mass trend 
% Last modified by bgetraer@princeton.edu, 1/5/2017
%**************************************************************************
addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
% Get 'thedates','ESTtotal','ESTtotalresid','total','alphavarall' from GREENLAND60.m
load(fullfile(datadir,'Greenland60data'));

%date domains:
Harig2013 = 1:monthnum(6,2013,thedates);
Getraer2018 = monthnum(6,2013,thedates):length(thedates);
before = 1:monthnum(6,2012,thedates);
during = monthnum(7,2012,thedates):monthnum(8,2013,thedates);
after = monthnum(9,2013,thedates):length(thedates);
all = 1:length(thedates);
without = [before after];
%% MODELS
%   y(t) = m(1) + m(2)*x + 1/2*m(3)*x^2 

% HARIG DATES
t = (thedates-thedates(1))/365; % time in years
s = total(Harig2013);
% 1st degree unweighted polynomial fit
[m1,f1,~,mint1] = linear_m(t,s,1);
% 2nd degree unweighted polynomial fit
[m2,f2,~,mint2] = linear_m(t,s,2);

% ALL DATES
tall = ((thedates-thedates(1))/365)'; % time in years
sall = total;

% 1st degree unweighted polynomial fit
[m1all,f1all,d1all,mint1all] = linear_m(tall,sall,1);

% 2nd degree unweighted polynomial fit
[m2all,f2all] = linear_m(tall,sall,2);

% AFTER DATES
tafter = ((thedates(after)-thedates(after(1)))/365)'; % time in years
safter = total(after);
s2after = total(after);

% 1st degree unweighted polynomial fit
[m1after,f1after] = linear_m(tafter,safter,1);

% 2nd degree unweighted polynomial fit
[m2after,f2after] = linear_m(tafter,safter,2);
[m21after,f21after] = linear_m(tafter,s2after,2);



%AGU 2016 polynomial
% acceleration2016=-28;

%% PLOTTING
%get positions for two subplots
figure(10)
for i = 1:2
    subplot(1,2,i)
    h = gca;
    pos(i,:) = h.Position;
end
close(gcf)

%set up some definitions
xlimit = [thedates(1)-100,thedates(end)+100];
ylimit = [-2500 1700];
yticklabels = -7:1:4;
yticksright = mmSLE2Gt(yticklabels);

ltextmodel_1 = strcat('\begin{tabular}{l}',...
    '1$^{st}$ order polynomial regression',...
    ' \end{tabular}'); %strcat('\def\du#1{\underline{\underline{#1}}}',...
%     '\begin{tabular}{l}',...
%     '1$^{st}$ order polynomial regression \\',... 
%     '$\underline{m}_{1}=(\du{G}^{T}_{1}*\du{G}_{1})^{-1}\du{G}^{T}_{1}*\underline{d}$',...
%     ' \end{tabular}');
ltextmodel_2 = strcat('\begin{tabular}{l}',...
    '2$^{nd}$ order polynomial regression',...
    ' \end{tabular}'); %strcat('\def\du#1{\underline{\underline{#1}}}',...
%     '\begin{tabular}{l}',...
%     '2$^{nd}$ order polynomial regression \\',... 
%     '$\underline{m}_{2}=(\du{G}^{T}_{2}*\du{G}_{2})^{-1}\du{G}^{T}_{2}*\underline{d}$',...
%     ' \end{tabular}');
ltextdata = strcat('\begin{tabular}{l}',...
    'monthly GRACE observations',...
    ' \end{tabular}');

% NOW START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure(1);
clf

% PLOT HARIG DATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax{1} = axes('Parent',f,'Position',pos(1,:));

% analysis of variance
resid1 = total(Harig2013)-f1(Harig2013);
resid2 = total(Harig2013)-f2(Harig2013);
vartest2(resid1,resid2)

hold on
errorbar(thedates,total,repmat(sqrt(alphavarall),size(total)),'color',[0.5 0.5 0.5],'linewidth',0.25);
totalplot = plot(thedates,total,'linewidth',1.5);
getraerplot = plot(thedates(Getraer2018),total(Getraer2018),'linewidth',1.5);
model1 = plot(thedates,f1,'k:','linewidth',2);
model2 = plot(thedates,f2,'--','linewidth',2);

%axes format
datetick
xlim(xlimit); ylim(ylimit);
xlabel('Year','interpreter','latex');set(gca,'fontsize',12);
ylabel('Mass (Gt)','interpreter','latex');

yyaxis right
ylabel('Sea Level Equivalence (mm)','interpreter','latex');

set(gca,'ylim',[yticksleft(1) yticksleft(end)],'ytick',yticksright,...
    'yticklabels',yticklabels,'ycolor',[0 0 0])

% accessory plotting
%new axes for text plotting
a = axes('Parent',f,'Position',pos(1,:));
axis off

%legend
lgd = legend([totalplot model1 model2] ,ltextdata,ltextmodel_1,ltextmodel_2);
set(lgd,'Interpreter','latex')


% text block of some interesting values
range = sprintf('Range $\\approx%i$ Gt',round(max(s)-min(s),-2));
m1slope = sprintf('$\\underline{m}_{1}$ Slope $=%0.2f\\pm{%0.2f}$ Gt per year',m1(2),mint1(2));
% Average modeled mass loss per year 2003-2013
    jan2003 = monthnum(1,2003,thedates);
    jan2013 = monthnum(1,2013,thedates);
    totalloss = f2(jan2013,:)-f2(jan2003,:);
modeledslope = sprintf('Avg. mass change 1/2003--1/2013 $=%i$ Gt per year',...
    round(totalloss/10));
m2acceleration = sprintf('$\\underline{m}_{2}$ Acceleration $=%0.2f\\pm{%0.2f}$ Gt per year$^{2}$',m2(3),mint2(3));
text(a,0.025,0.1,...
    sprintf('\\begin{tabular}{l} %s %s %s %s %s %s %s \\end{tabular}',...
    range,'\\',m1slope,'\\',modeledslope,'\\',m2acceleration),'interpreter','latex','fontsize',12)

% title
line1 = 'Greenland GRACE signal, 2003--2014';
line2 = 'Recreation of Harig~\&~Simons, 2016, their Figure 4c';

title(sprintf('\\begin{tabular}{c} \\textbf{%s} %s %s \\end{tabular}',line1,'\\',line2),...
    'interpreter','latex','fontsize',12,'horizontalalignment','center')


% ALL DATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax{2} = axes('Parent',f,'Position',pos(2,:));

hold on
errorbar(thedates,total,repmat(sqrt(alphavarall),size(total)),'color',[0.5 0.5 0.5],'linewidth',0.25);
total_m = plot(thedates,total,'linewidth',1.5);
% total_e = plot(thedates,ESTtotal,'linewidth',1.5);
model1 = plot(thedates,f1all,'k:','linewidth',2);
         plot(thedates,f1all-d1all,'k:','linewidth',0.5);
delta1 = plot(thedates,f1all+d1all,'k:','linewidth',0.5);
% model2 = plot(thedates,f2all,'--','linewidth',2);
% analysis of variance
resid1 = total(:)-f1all(:);
resid2 = total(:)-f2all(:);
vartest2(resid1,resid2)
[r p] = corrcoef(total,f1all)

%axes format
datetick
xlim(xlimit); ylim(ylimit);
xlabel('Year','interpreter','latex');set(gca,'fontsize',12);
ylabel('Mass (Gt)','interpreter','latex');

yyaxis right
ylabel('Sea Level Equivalence (mm)','interpreter','latex');

set(gca,'ylim',[yticksleft(1) yticksleft(end)],'ytick',yticksright,...
    'yticklabels',yticklabels,'ycolor',[0 0 0])

% accessory plotting
%new axes for text plotting
a = axes('Parent',f,'Position',pos(2,:));
axis off

% text block of some interesting values
range = sprintf('Range $\\approx%i$ Gt',round(max(sall)-min(sall),-2));
m1slope = sprintf('$\\underline{m}_{1}$ Slope $=%0.2f\\pm{%0.2f}$ Gt per year',m1all(2),mint1all(2));
% Average modeled mass loss per year 2003-2013
    jan2003 = monthnum(1,2003,thedates);
    jan2016 = monthnum(1,2013,thedates);
    totalloss = f1all(jan2016,:)-f1all(jan2003,:);
modeledslope = sprintf('Avg. mass change 1/2003--1/2017 $=%i$ Gt per year',...
    round(totalloss/10));
rsq = sprintf('Correlation of model: R$^2=%0.3f$',r(2));
text(a,0.025,0.1,...
    sprintf('\\begin{tabular}{l} %s %s %s %s %s %s %s \\end{tabular}',...
    range,'\\',m1slope,'\\',modeledslope,'\\',rsq),'interpreter','latex','fontsize',12)

%legend
lgd = legend([total_m model1 ] ,ltextdata,ltextmodel_1);
set(lgd,'Interpreter','latex')

% title
line1 = 'Greenland GRACE signal, 2003--2017';
title(sprintf('\\begin{tabular}{c} \\textbf{%s} \\end{tabular}',line1),...
    'interpreter','latex','fontsize',12,'horizontalalignment','center')