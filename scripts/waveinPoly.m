%**************************************************************************
%   Script developed for "REGIONAL FORCING OF GREENLAND ICE LOSS 2002-2017"
%   Spring 2018 Junior Paper, Princeton Department of Geosciences
%
%   SCRIPT 5
%   
%   PREVIOUS: ANALYZEWAVELET.m
%   NEXT: .m
%
%   Benjamin Getraer bgetraer@princeton.edu
%   Modified: 2/22/2019
%   See also: 
%**************************************************************************

% locate slepian_bgetraer function and datafile directories, and set workspace
homedir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/';
functiondir = fullfile(homedir,'functions');
datadir = fullfile(homedir,'datafiles');
addpath(functiondir,datadir);   clear('homedir','functiondir');
setworkspace();

% load datafiles created in BOXGREENLAND.m and IMAGERYSEQ.m
load ptsGL 
load im_tools 
load im_seqSH
%% Spatial index for points in and around Greenland

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MASK OF GREENLAND IN IMAGE BASIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz = size(xp);
a = inpolygon(xp(:),yp(:),bx,by);   % array of points inside Greenland
A = 1*reshape(a,sz);                % matrix of points inside Greenland

% plot a visualization
figure(1)
clf
imshow(A)
hold on
title("mask of image points inside of Greenland")
colormap(bone)
plot(bx,by,'r','linewidth',2)
plot(gx,gy,'k','linewidth',1)

save(fullfile(datadir,'GLimagemask'),'A')
%% Visualization of level indexing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   IMAGE RECONSTRUCTION BY DECOMPOSITION LEVEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
wavelevelINDEX('test');

%% Threshold by spatial support within Greenland

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   GENERATE UNIFORM WAVELET DECOMPOSITION DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavename = 'haar';
level = wmaxlev(sz,wavename);
[C,S]=wavedec2(A,level,wavename);
C0 = zeros(1,length(C));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FIND THE WAVELETS WHICH PASS THE THRESHOLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ I, L, coefinlevel] = wavelevelINDEX( C,S );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the area threshold a wavelet of a given level must occupy within the
% polygon (by percent of wavelet support) in order to pass.
%   highest resolution, 90%
%   top two lowest resolution, take them all
areathresh = [linspace(0.9,0,level-1) 0];

% initialize pass array 
polygonpassindex = zeros(size(C));
if lengthyProcessFlag('wavelet in polygon index')
    for i = 1:level
        fprintf('level %d \n',i)
        % if we are keeping them all anyways!
        if areathresh(i)==0
            polygonpassindex(I{i})=1;
        else
            % only look at one of the wavelet sequences
            index = I{i}(1:coefinlevel(i));
            % passes per level
            passlevel = zeros(size(index));
            for in = index
                fprintf('wavelet %d of %d \n',find(index==in),coefinlevel(i))
                Cmod = C0;
                Cmod(in) = 1;
                imj = waverec2(Cmod,S,wavename);
                xw = xp(imj(:)~=0);
                yw = yp(imj(:)~=0);
                pw = inpolygon(xw(:),yw(:),bx,by);
                aw = sum(imj(:)~=0);
                ainp = sum(pw(:)~=0);
                passlevel(index==in) = (1-(aw-ainp)/aw > areathresh(i));
            end
            
            polygonpassindex(I{i}) = [passlevel, passlevel, passlevel];
        end
    end
    
    filename = 'inGLpass';
    save(fullfile(datadir,filename),'polygonpassindex')
end
%%
figure(6)
clf
% [wd,s]=wavedec2(A,level,wavename);
PASS = waverec2(polygonpassindex,S,wavename);
imagesc(PASS)
colormap(bone);
axis image
%%
filename = sprintf('pass%dlin2',order);
load(fullfile(datadir,filename))
figure(6)
clf
% [wd,s]=wavedec2(A,level,wavename);
PASS = waverec2(polygonpassindex,s,wavename);
imagesc(PASS)
colormap(bone);
axis image
hold on;
% Greenland
plot(bx,by,'w-')
axis off
pause(0.01)

%% HERE IS A MONEY FIGURE
Ddiff = D(:,:,end) - D(:,:,1);
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
wDTA = waverec2(DT.*polygonpassindex,sdiff,wavename);
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
wDTA = waverec2(DT.*polygonpassindex,sdiff,wavename);
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
AT = waverec2(wA.*polygonpassindex,sdiff,wavename);
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

%%

bias3 = sprintf('\\textbf{bias} & = &%0.3e',abs((mean(Ddiff(:).*A(:))-mean(wDTA(:).*A(:)))/mean(Ddiff(:))))
bias3 = sprintf('\\textbf{bias} & = &%0.3e',abs((sum(Ddiff(:).*A(:))-sum(wDTA(:).*A(:)))/sum(Ddiff(:))))

%%
% invariance curve for the classified image
level=8;
wavename = 'haar';
[c,s]=wavedec2(A,level,wavename);

pthresh = 0.99

figure(2)
clf
hold on
c = 'r';
[wdiff,sA]=wavedec2(Ddiff,level,wavename);
abwA = abs(wdiff);
N = 1:level;
ptl = linspace(90,100,50);
Ddiff =  D(:,:,end)-D(:,:,1);

clear T er r2 wD b NC b
for i = 1:length(ptl)
    T(i) = prctile(abwA,ptl(i));
    NC{i} = wthcoef2('t',wdiff,sA,N,repmat(T(i),size(N)),'h');
    imj{i} = waverec2(NC{i},sA,wavename);
    er(i) = immse(Ddiff,imj{i});
    r2(i) = 1-var(Ddiff(:)-imj{i}(:))/var(Ddiff(:));
    b(i) =  abs((sum(Ddiff(:))-sum(imj{i}(:)))/sum(Ddiff(:)));
end
hb = plot(ptl,b,':','color',c,'linewidth',2);
h = plot(ptl,r2,'color',c,'linewidth',2);
[~, qc] = min(abs(r2-pthresh));  % where does the curve reach 90% invariance
QC = wthcoef2('t',wdiff,sA,N,repmat(T(qc),size(N)),'h');
np = sum(QC~=0);
Q = waverec2(QC,sA,wavename);

txth=strcat('\begin{tabular}{lr} \textbf{',wavename,'} inv. &',...
    num2str(np), sprintf(' coefficients at %0.2f',pthresh), '\end{tabular}');
txtb=strcat('\begin{tabular}{lr} \textbf{',wavename,'} bias &',...
    sprintf('%0.3e  at %0.2f',b(qc),pthresh), '\end{tabular}');

axis([90 100 ylim])
plot(xlim,[pthresh pthresh],'--','color',0.6*[1 1 1],'linewidth',1)
ylabel('image invariance or bias')
xlabel('percentile threshold')
lgd = legend([h, hb],txth, txtb);
set(lgd,'interpreter','latex')
grid on
%% CHOOSE HAAR WAVELET, CHOOSE percentile threshold to minimize bias
wavename = 'haar';
ptile = 99.6;%ptl(qc);
level = 8;

[wA,sA]=wavedec2(A,level,wavename);
abwA = abs(wA);
N = 1:level;
clear T er r2 wD b NC b
T = prctile(abwA,ptile);
wAT = wthcoef2('t',wA,sA,N,repmat(T,size(N)),'h');
AT = waverec2(wAT,sA,wavename);
r2 = 1-var(Ddiff(:)-AT(:))/var(Ddiff(:));
b =  abs((sum(Ddiff(:))-sum(AT(:)))/sum(Ddiff(:)));
%
% hb{j} = plot(ptl,b,':','color',c{j},'linewidth',2);
% h{j} = plot(ptl,r2,'color',c{j},'linewidth',2);
% [~, qc(j)] = min(abs(r2-pthresh));  % where does the curve reach 90% invariance
% QC = wthcoef2('t',wdiff,sdiff,N,repmat(T(qc(j)),size(N)),'h');
% np(j) = sum(QC~=0);
% Q{j} = waverec2(QC,sdiff,wname{j});


figure(3)
clf
wcoff = sort(abwA,'descend');
X = log(wcoff);
% hg = histogram(log(wcoff),200);

hold on

%// Histogram plot:
[y n] = hist(X(X>-60),200); %// y: values; n: bin centers
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
sum(wAT~=0)


%%
% a plot of how well these wavelets capture greenland
figure(3)
subplot(2,1,1)
imagesc(A)
axis image
hold on
plot(gx,gy,'k-')
colormap(hot);
caxis([0 1])

subplot(2,1,2)
imagesc(AT)
axis image
hold on
plot(gx,gy,'k-')
caxis([0 1])
colorbar

save(fullfile(datadir,'wavemask'),'wAT')
