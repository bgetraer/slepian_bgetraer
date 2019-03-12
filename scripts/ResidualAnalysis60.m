%**************************************************************************
% Integrated Mass trend 
% Last modified by bgetraer@princeton.edu, 1/3/2017
%**************************************************************************
addpath('/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/functions')
setworkspace('/Users/benjamingetraer/Documents/IndependentWork/SH_Workspace');
datadir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/datafiles';
% Get 'thedates','ESTtotal','ESTtotalresid','total','alphavarall' from GREENLAND60.m
load(fullfile(datadir,'Greenland60data'));

%date domains:
Harig2013 = 1:monthnum(6,2013,thedates);
Getraer2018 = monthnum(6,2013,thedates):length(thedates);
before = 1:monthnum(6,2012,thedates);
during = monthnum(7,2012,thedates):monthnum(8,2013,thedates);
after = monthnum(9,2013,thedates):length(thedates);
all = 1:length(thedates);

% errorbar values
twosigma = 2*sqrt(alphavarall);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMINE RESIDUALS
f = figure(3);
clf
subplot(1,2,1)
hold on
ebarhdl = errorbar(thedates,total,repmat(twosigma,size(total)),'color',[0.5 0.5 0.5],'linewidth',0.25);
s = ebarhdl.LineStyle;ebarhdl.LineStyle = 'none';
ESTtotalhdl = plot(thedates,ESTtotal,'linewidth',1.5);
totalhdl = plot(thedates,total,'linewidth',1.5);

%set up some definitions
xlimit = [thedates(1)-100,thedates(end)+100];
ylimit = [-2000 1800];
yticklabels = -7:1:4;
yticksright = mmSLE2Gt(yticklabels);

datetick; grid on;
xlim(xlimit); ylim(ylimit);
xlabel('Year','interpreter','latex');set(gca,'fontsize',12);
ylabel('Mass (Gt)','interpreter','latex');
yticksleft = get(gca,'ytick');

yyaxis right
ylabel('Sea Level Equivalence (mm)','interpreter','latex');

set(gca,'ylim',[yticksleft(1) yticksleft(end)],'ytick',yticksright,...
    'yticklabels',yticklabels,'ycolor',[0 0 0])
%legend
lgd = legend([totalhdl,ESTtotalhdl,ebarhdl],'GRACE signal','$f_{\alpha}$ model','$2\sigma$');
set(lgd,'Interpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subax = axes('Parent',f,'Position',[0.1300 0.1100 0.15  0.4]);
axis off
buff=greenland(10,0.5);
gl=greenland(10,0);
XY1 = [buff;gl;];
XY2 = gl;
[x1,y1,z1] = sph2cart(XY1(:,1)*pi/180,XY1(:,2)*pi/180,1); 
[x2,y2,z2] = sph2cart(XY2(:,1)*pi/180,XY2(:,2)*pi/180,1); 
hold on
patch(x2,y2,z2,[1,1,1]); patch(x1,y1,z1,[1,0,0.1]);
view(45,90)
axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2)
hold on
resid_e = plot(thedates,ESTtotalresid,'-o','linewidth',1.5,'markerfacecolor',...
    [0.7 0 0.1],'markeredgecolor','none');
% standards of deviation
[sigma_all,~,sigma_allhdl] = stdline(thedates,ESTtotalresid);

% defined as before the 2012 excursion below sigma_all
% sigma_before = stdline(thedates(before),ESTtotalresid(before),'r','-');
% sigma_after = stdline(thedates(after),ESTtotalresid(after),'r','-');
% sigma_during = stdline(thedates(during),ESTtotalresid(during),'r','-');

%axes format
datetick; grid on;
xlim([thedates(1)-100,thedates(end)+100]); ylim([-300 300]);
xlabel('Year');ylabel('Mass (Gt)');set(gca,'fontsize',12);

yticklabels = -1.5:0.25:1.5;
yticksright = mmSLE2Gt(yticklabels);

xlabel('Year','interpreter','latex');set(gca,'fontsize',12);
ylabel('Mass (Gt)','interpreter','latex');
yticksleft = get(gca,'ytick');

yyaxis right
ylabel('Sea Level Equivalence (mm)','interpreter','latex');

set(gca,'ylim',[yticksleft(1) yticksleft(end)],'ytick',yticksright,...
    'yticklabels',yticklabels,'ycolor',[0 0 0])

% plot the 2std errorbar values
plot(xlim,[twosigma twosigma],'-','color',[0.25 0.25 0.25],'linewidth',1);
twosigma_allhdl = plot(xlim,[-twosigma -twosigma],'-','color',[0.25 0.25 0.25],'linewidth',1);

%legend
lgd = legend([resid_e,sigma_allhdl,twosigma_allhdl],'residuals','$\sigma$','$2\sigma$');
set(lgd,'Interpreter','latex')

%%
% y=ESTtotalresid(before)/sigma_before;
% y=ESTtotalresid(after)/sigma_after;
% y=ESTtotalresid(all)/sigma_all;
y=ESTtotalresid(during)/sigma_during;


y=y-mean(y);
x = round(min(y)-1):0.1:round(max(y)+1);
figure(1)
clf
h=histogram(y,'binmethod','scott','normalization','pdf');
hold on
plot(x,normpdf(x,0,1));
[H P] =kstest(y)

clf
[f,x_values] = ecdf(y);
F = plot(x_values,f);
set(F,'LineWidth',2);
hold on;
G = plot(x_values,normcdf(x_values,0,1),'r-');