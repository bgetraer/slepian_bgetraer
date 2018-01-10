%% SETUP WORKSPACE AND IMPORT SPHERICAL HARMONIC DATA
addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

% import all data in the GRACE original directory and organize as plmt 
% potcoffs: matrix of SH coefficients for the geopotential field 
% cal_errors: estimated errors given with the data
%   format: (monthnumber,l_m,values)
% thedates: vector of datenum associated with monthnumber
%   thedates are midpoints of data collection times

[potcoffs,cal_errors,thedates]=grace2plmt('CSR','RL05','SD',0);
%% May 2003
figure(5)
plotplm(squeeze(potcoffs(monthnum(5,2003,thedates),:,:)),[],[],[],0.1);
caxis([-2 2]*1e6)
c1 = colorbar;
c1.Label.String='kg per m$^2$';
c1.Label.Interpreter = 'latex';
title('\textbf{Surface mass density derived from Geopotential field, May 2003}','Interpreter','latex')
%% Compare October and April
% difference coeff matrix
year = 2007;

theelsems=squeeze(potcoffs(1,:,1:2));   % preserve the lm section
Oct = squeeze(potcoffs(monthnum(10,year,thedates),:,3:4));
Apr = squeeze(potcoffs(monthnum(4,year+1,thedates),:,3:4));
diff_coef = [theelsems Apr-Oct];

figure(3)
plotplm(diff_coef,[],[],2,0.2)
caxis([-.6 1]*1e3)
c2 = colorbar;
c2.Label.String='surface mass density';
title('Difference of Apr. 2008 and Oct. 2007, unfiltered')

 figure(3);hold on;
filter_L = 30;
plotplm(plmfilt(diff_coef,filter_L))
caxis([-.4 .8]*1e3)
c3 = colorbar;
c3.Label.String='surface mass density';
title(sprintf('Difference of Apr. 2008 and Oct. 2007, plmfilt(L=%i)',filter_L))
%% Compare February and August
% difference coeff matrix
year = 2008;

theelsems=squeeze(potcoffs(1,:,1:2));   % preserve the lm section
Aug = squeeze(potcoffs(monthnum(8,year,thedates),:,3:4));
Feb = squeeze(potcoffs(monthnum(2,year,thedates),:,3:4));
diff_coef = [theelsems Aug-Feb];

f =figure(1)
subplot(1,2,1)
plotplm(diff_coef,[],[],5)
caxis([-.6 1]*1e3)
title('Unfiltered','interpreter','latex')
set(gca,'fontsize',12)


subplot(1,2,2);
hold on;
filter_L = 30;
plotplm(plmfilt(diff_coef,filter_L),[],[],5)
caxis([-.6 1]*1e3)
title('Filtered','interpreter','latex')
set(gca,'fontsize',12)


ax = copyobj(gca,f)
axis off
c2 = colorbar;
c2.Label.String = 'kg per m$^2$';
c2.Position = [0.91 0.2667 0.0299 0.5019];
c2.Label.Interpreter = 'latex';
set(gca,'fontsize',12)

%% TEST WITH plotplm.m FROM slepian_alpha

figure(1)
for index=2
% plot all orders and degrees from MAY 2008
theelsems=squeeze(potcoffs(1,:,1:2));
plotplm([theelsems squeeze(potcoffs(index,:,3:4))-squeeze(potcoffs(index-1,:,3:4))])
title(sprintf('month %i',index))
pause
end
% figure(2)
% plot order 1 from MAY 2008
% plotplm(squeeze(potcoffs(1,2,:))')

%% show second order terms
% Create blank LMCOSI matrix
[~,~,~,blank,~,~,~,~,~,ronm]=addmon(2);
blank(blank(:,1)==2&blank(:,2)==1,3:4) = [1,1];
% blank(blank(:,1)==2&blank(:,2)==0,3:4) = [1,0];
blank

plotplm(blank,[],[],1,1)