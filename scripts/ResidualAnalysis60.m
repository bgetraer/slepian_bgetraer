%**************************************************************************
% Integrated Mass trend 
% Last modified by bgetraer@princeton.edu, 1/3/2017
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
%%
%% EXAMINE RESIDUALS
figure(3)
clf
hold on
resid_e = plot(thedates,ESTtotalresid,'linewidth',1.5);
% standards of deviation
sigma_all = stdline(thedates,ESTtotalresid);
% defined as before the 2012 excursion below sigma_all
sigma_before = stdline(thedates(before),ESTtotalresid(before),'r','-');
sigma_after = stdline(thedates(after),ESTtotalresid(after),'r','-');
sigma_during = stdline(thedates(during),ESTtotalresid(during),'r','-');

%axes format
datetick
xlim([thedates(1)-100,thedates(end)+100]); ylim([-500 500]);
xlabel('Year');ylabel('Mass (Gt)');set(gca,'fontsize',12);
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