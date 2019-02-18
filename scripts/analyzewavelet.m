%**************************************************************************
%   Script developed for "REGIONAL FORCING OF GREENLAND ICE LOSS 2002-2017"
%   Spring 2018 Junior Paper, Princeton Department of Geosciences
%
%   SCRIPT 4
%   
%   PREVIOUS: CHOOSEWAVELET.m
%   NEXT: WAVEINPOLY.m
%
%   Benjamin Getraer bgetraer@princeton.edu
%   Modified: 7/25/2018
%   See also: 
%**************************************************************************

% locate slepian_bgetraer function and datafile directories, and set workspace
homedir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/';
functiondir = fullfile(homedir,'functions');
datadir = fullfile(homedir,'datafiles');
addpath(functiondir,datadir);   clear('homedir','functiondir');
setworkspace();


% load datafiles created in BOXGREENLAND.m and IMAGERYSEQ.m
load('ptsGL')
load('im_tools')
load('im_seqSH')

order = 10;
%% PERFORM DECONSTRUCTION ON THE FULL TIMESERIES DIFFERENCE IMAGE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Caluculate coefficient percentile threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the difference between the first and last date in the timeseries
Ddiff =  D(:,:,end)-D(:,:,1);
wavename = 'haar';
invarT = 0.90; % invariance threshold
% this is the percentile at 90% invar
[ptile, level] = prctileThold(Ddiff, invarT,wavename); 




[wdiff,sdiff]=wavedec2(Ddiff,level,wavename);
abwdiff = abs(wdiff);
N = 1:level;
        clear T er r2 wD b NC b
T = prctile(abwdiff,ptile);
DT = wthcoef2('t',wdiff,sdiff,N,repmat(T,size(N)),'h');
wDT = waverec2(DT,sdiff,wavename);
r2 = 1-var(Ddiff(:)-wDT(:))/var(Ddiff(:));
b =  abs((sum(Ddiff(:))-sum(wDT(:)))/sum(Ddiff(:)));
% 
% hb{j} = plot(ptl,b,':','color',c{j},'linewidth',2);
% h{j} = plot(ptl,r2,'color',c{j},'linewidth',2);
% [~, qc(j)] = min(abs(r2-pthresh));  % where does the curve reach 90% invariance
% QC = wthcoef2('t',wdiff,sdiff,N,repmat(T(qc(j)),size(N)),'h');
% np(j) = sum(QC~=0);
% Q{j} = waverec2(QC,sdiff,wname{j});
    
    
figure(2)
clf
wcoff = sort(abwdiff,'descend');
X = log(wcoff);
% hg = histogram(log(wcoff),200);

hold on
% set(gca, 'xticklabels', exp(xticks))
% plot(log(wcoff))

%// Histogram plot:
[y n] = hist(X,200); %// y: values; n: bin centers
yn = y/sum(y);
lT = log(T);
ind = n>lT; %// bin centers: greater or smaller than threshold?
bar(n(ind), yn(ind), 1, 'b'); %// for greater: use red
hold on %// keep graph, Or use hold(your_axis_handle, 'on')
bar(n(~ind), yn(~ind), 1, 'r','edgecolor','k'); %// for smaller: use blue

[~, nd] = min(abs(n-lT)); %// locate bar around D: it needs the two colors
patch([(n(nd+1)+n(nd))/2 lT lT (n(nd+1)+n(nd))/2], [0 0 yn(nd) yn(nd)], 'b');
%// take care of that bar with a suitable patch
plot(log([T T]),ylim,'k')

% percentage of coefficients thrown out:
sum(yn(~ind))
% number left
sum(DT(:)~=0)

threshpass = DT~=0;
filename = sprintf('thresh%d',order);
save(fullfile(datadir,filename),'threshpass')
%% MOVED TO WAVEINPOLY
figure(3)
clf
subplot(3,1,1)
imagesc(Ddiff)
axis image
hold on
plot(gx,gy,'k-')
colormap(bluewhitered(1000,1));
title(sprintf('%d wavelets',sum(wdiff~=0)))
text1 = sprintf('bias = %0.3e',abs((sum(Ddiff(:))-sum(Ddiff(:)))/sum(Ddiff(:))));
text(300,150,text1)
axis off

subplot(3,1,2)
imagesc(wDT)
axis image
hold on
plot(gx,gy,'k-')
colormap(bluewhitered(1000,1));
title(sprintf('%d wavelets',sum(DT~=0)))
text2 = sprintf('bias = %0.3e',abs((sum(Ddiff(:))-sum(wDT(:)))/sum(Ddiff(:))));
text(300,150,text2)
axis off

subplot(3,1,3)
wDTA = waverec2(DT.*pass,sdiff,wavename);
imagesc(wDTA)
axis image
hold on
plot(gx,gy,'k-')
colormap(bluewhitered(1000,1));
title(sprintf('%d wavelets',sum(DT.*pass~=0)))
bias3 = sprintf('bias = %0.3e',abs((sum(Ddiff(:))-sum(wDTA(:)))/sum(Ddiff(:))));
invar3 = sprintf('invar = %0.3f',1-var(Ddiff(:)-wDTA(:))/var(Ddiff(:)));
text3 = strcat('\begin{tabular}{l}',invar3,'\\',bias3,'\end{tabular}');
text(300,150,text3,'interpreter','latex')
axis off

% % f = 1:9
% % F = reshape(f,[3 3])
% T = prctile(sortwaveC,90)
% N = [1:level];
% NC = wthcoef2('t',wdiff,sdiff,N,repmat(T,size(N)),'h');
% 
% figure(4)
% wDT = waverec2(NC,sdiff,wavename);
% imagesc(wDT)
% axis image
% hold on
% plot(gx,gy,'k-')
% colormap(bluewhitered(1000,1));