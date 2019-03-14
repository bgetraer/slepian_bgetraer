function jp02fig7()
%JP02FIG7 Creates Figure 7 from "REGIONAL FORCING OF GREENLAND ICE LOSS 
%   2002-2017" Spring 2018 Junior Paper, Princeton Department of Geosciences
%
% last modified by: bgetraer@princeton.edu, 3/12/2019

homedir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/';
functiondir = fullfile(homedir,'functions');
datadir = fullfile(homedir,'datafiles');
addpath(functiondir,datadir);   clear('homedir','functiondir');
load ptsGL 
load im_tools 
load im_seqSH
load threshpassindex

wavename = 'haar';
Ddiff = D(:,:,end) - D(:,:,1);
sz = size(Ddiff);
level = wmaxlev(sz,wavename);
[wdiff,S]=wavedec2(Ddiff,level,wavename);


filename = 'passIndexbyAreaONLY.mat';
if ~exist(fullfile(datadir,filename),'file')
    [ polygonpassindexONLY ] = ...
        areaInPolygon( wdiff, S, wavename, level, xp, yp, bx, by );
    filename = 'passIndexbyAreaONLY';
    save(fullfile(datadir,filename),'polygonpassindexONLY')
end
load(fullfile(datadir,filename))
filename = 'passIndexbyArea';
load(fullfile(datadir,filename))
load(fullfile(datadir,'GLimagemask'))


 
DT = wdiff.*threshpassindex;
wDT = waverec2(DT,S,wavename);

figure(3)
clf
subplot(3,3,1)
imagesc(Ddiff)
axis image
hold on
plot(gx,gy,'k-')
title(sprintf('%d wavelets',sum(wdiff~=0)))
bias1 = sprintf('\\textbf{bias} & = &%0.3e',abs((sum(Ddiff(:))-sum(Ddiff(:)))/sum(Ddiff(:))));
invar1 = sprintf('\\textbf{invar} &= &%0.3f',1-var(Ddiff(:)-Ddiff(:))/var(Ddiff(:)));
text1 = strcat('\begin{tabular}{lcr}',invar1,'\\',bias1,'\end{tabular}');
% text(300,150,text1,'interpreter','latex')
axis off
cax = caxis;


subplot(3,3,4)
imagesc(wDT)
axis image
hold on
plot(gx,gy,'k-')
title(sprintf('%d wavelets',sum(DT~=0)))
bias2 = sprintf('\\textbf{bias} & = &%0.1e',abs((sum(Ddiff(:))-sum(wDT(:)))/sum(Ddiff(:))));
invar2 = sprintf('\\textbf{invar} &= &%0.3f',1-var(Ddiff(:)-wDT(:))/var(Ddiff(:)));
text2 = strcat('\begin{tabular}{lcr}',invar2,'\\',bias2,'\end{tabular}');
text(-80,128,text2,'interpreter','latex','rotation',90,'horizontalalignment','center')
caxis(cax)
axis off

subplot(3,3,7)
wDTA = waverec2(DT.*polygonpassindex,S,wavename);
imagesc(wDTA)
axis image
hold on
plot(gx,gy,'k-')
% text and titles
title(sprintf('%d wavelets',sum(DT.*polygonpassindex~=0)))
bias3 = sprintf('\\textbf{bias} & = &%0.1e',abs((sum(Ddiff(:))-sum(wDTA(:)))/sum(Ddiff(:))));
invar3 = sprintf('\\textbf{invar} &= &%0.3f',1-var(Ddiff(:)-wDTA(:))/var(Ddiff(:)));
text3 = strcat('\begin{tabular}{lcr}',invar3,'\\',bias3,'\end{tabular}');
text(-80,128,text3,'interpreter','latex','rotation',90,'horizontalalignment','center')
axis off
caxis(cax)

% NOW FOR THE MASKED IMAGES
subplot(3,3,3)
imagesc(Ddiff.*A)
axis image
hold on
plot(gx,gy,'k-')
title(sprintf('%d wavelets, masked',sum(wdiff~=0)))
axis off
caxis(cax)

subplot(3,3,6)
imagesc(wDT.*A)
axis image
hold on
plot(gx,gy,'k-')
title(sprintf('%d wavelets, masked',sum(DT~=0)))
bias2 = sprintf('\\textbf{bias} & = &%0.2f',abs((sum(Ddiff(:).*A(:))-sum(wDT(:).*A(:)))/sum(Ddiff(:).*A(:))));
invar2 = sprintf('\\textbf{invar} &= &%0.3f',1-var(Ddiff(:).*A(:)-wDT(:).*A(:))/var(Ddiff(:).*A(:)));
text2 = strcat('\begin{tabular}{lcr}',invar2,'\\',bias2,'\end{tabular}');
text(256+80,128,text2,'interpreter','latex','rotation',90,'horizontalalignment','center')
axis off
caxis(cax)

subplot(3,3,9)
wDTA = waverec2(DT.*polygonpassindex,S,wavename);
imagesc(wDTA.*A)
axis image
hold on
plot(gx,gy,'k-')
% text and titles
title(sprintf('%d wavelets, masked',sum(DT.*polygonpassindex~=0)))
bias3 = sprintf('\\textbf{bias} & = &%0.3f',abs((sum(Ddiff(:).*A(:))-sum(wDTA(:).*A(:)))/sum(Ddiff(:).*A(:))));
invar3 = sprintf('\\textbf{invar} &= &%0.3f',1-var(Ddiff(:).*A(:)-wDTA(:).*A(:))/var(Ddiff(:).*A(:)));
text3 = strcat('\begin{tabular}{lcr}',invar3,'\\',bias3,'\end{tabular}');
text(256+80,128,text3,'interpreter','latex','rotation',90,'horizontalalignment','center')
axis off
caxis(cax)

% MIDDLE IMAGES

subplot(3,3,2)
imagesc(cax(1)*A)
axis image
hold on
plot(gx,gy,'w-')
% text and titles
title(sprintf('buffered wavelet mask',sum(DT.*polygonpassindex~=0)))
axis off
caxis(cax)

subplot(3,3,5)
[wA,s]=wavedec2(A,level,wavename);
AT = waverec2(wA.*polygonpassindexONLY,S,wavename);
imagesc(cax(1)*AT)
axis image
hold on
plot(gx,gy,'w-')
% text and titles
title(sprintf('threshold by area',sum(DT.*polygonpassindex~=0)))
bias4 = sprintf('\\textbf{bias} & = &%0.1e',abs((sum(A(:))-sum(AT(:)))/sum(A(:))));
invar4 = sprintf('\\textbf{invar} &= &%0.3f',1-var(A(:)-AT(:))/var(A(:)));
text4 = strcat('\begin{tabular}{lcr}',invar4,'\\',bias4,'\end{tabular}');
text(128,290,text4,'interpreter','latex','horizontalalignment','center')
axis off
caxis(cax)

subplot(3,3,8)

ttl = sprintf('%s--%s',datestr(thedates(1),'mmmm yyyy'),datestr(thedates(end),'mmmm yyyy'));

text(0.5,0.7,'\bfseries{Greenland Mass Wasting}','interpreter','latex',...
    'horizontalalignment','center','fontsize',14)
text(0.5,0.5,ttl,'interpreter','latex','horizontalalignment','center',...
    'fontsize',12)

axis off
caxis(cax)
cb = colorbar;
set(cb,'Location', 'south','AxisLocation','out','ticks',[-3500 -1750 0 1500])
ylabel(cb,'kg per m^2');

colormap(bluewhitered(1000,1));
end

