addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

order = 10;
load(strcat('ptsGL',num2str(order)))
load(strcat('im_tools',num2str(order)))
load(strcat('im_endsSH',num2str(order)))

close all
%%
% find the points inside of Greenland
a = inpolygon(xp(:),yp(:),bx,by);
A = 1*reshape(a,size(xp));

% plot
figure(1)
clf
imagesc(A)
hold on
colormap(bone)
plot(bx,by,'r','linewidth',2)
plot(gx,gy,'k','linewidth',1)
%% where are the different level?
level = 8;
levelin{8} = 1:4*l(8);
for i = (length(l)-1):-1:1
    % level i
    levelin{i} = levelin{i+1}(end) + (1:3*l(i));
end
%%
% DETERMINE WHICH WAVELETS PASS AN AREA THRESHOLD WITHIN A BUFFERED
% GREENLAND
[wd,s]=wavedec2(A,level,wavename);
wd = wd*0+1;

areathresh = linspace(0.9,0.05,level-1);

pass = zeros(size(wd));

for i = 1:level-1
    fprintf('level %d \n',i)
    
    index = levelin{i}(1:l(i)); % only look at one of the wavelet sequences
    % per level
    passlevel = zeros(size(index));
    for in = index
        %         fprintf('wavelet %d of %d \n',find(index==in),l(i))
        thiswd = wd;
        thiswd([1:in-1,in+1:end]) = 0;
        wD = waverec2(thiswd,s,wavename);
        xw = xp(wD(:)~=0);
        yw = yp(wD(:)~=0);
        pw = inpolygon(xw(:),yw(:),bx,by);
        aw = sum(wD(:)~=0);
        ainp = sum(pw(:)~=0);
        passlevel(index==in) = (1-(aw-ainp)/aw > areathresh(i));
%         % optional plotting
%         imagesc(~wD)
%         colormap(bone);
%         axis image
%         hold on;
%         % Greenland
%         plot(gx,gy,'k-')
%         axis off
%         pause(0.1)
        
    end
    
    pass(levelin{i}) = [passlevel, passlevel, passlevel];
end

pass(levelin{8})=1;
pass(levelin{7})=1;

filename = sprintf('pass%dlin',order);
save(fullfile(datadir,filename),'pass')
%%
load(fullfile(datadir,filename))
figure(6)
clf
% [wd,s]=wavedec2(A,level,wavename);
PASS = waverec2(pass,s,wavename);
imagesc(PASS)
colormap(bone);
axis image
hold on;
% Greenland
plot(gx,gy,'k-')
axis off
pause(0.01)

%% HERE IS A MONEY FIGURE

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
bias2 = sprintf('\\textbf{bias} & = &%0.3e',abs((sum(Ddiff(:))-sum(wDT(:)))/sum(Ddiff(:))));
invar2 = sprintf('\\textbf{invar} &= &%0.3f',1-var(Ddiff(:)-wDT(:))/var(Ddiff(:)));
text2 = strcat('\begin{tabular}{lcr}',invar2,'\\',bias2,'\end{tabular}');
text(-80,128,text2,'interpreter','latex','rotation',90,'horizontalalignment','center')
caxis(cax)
axis off

subplot(3,3,7)
wDTA = waverec2(DT.*pass,sdiff,wavename);
imagesc(wDTA)
axis image
hold on
plot(gx,gy,'k-')
% text and titles
title(sprintf('%d wavelets',sum(DT.*pass~=0)))
bias3 = sprintf('\\textbf{bias} & = &%0.3e',abs((sum(Ddiff(:))-sum(wDTA(:)))/sum(Ddiff(:))));
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
bias2 = sprintf('\\textbf{bias} & = &%0.3e',abs((sum(Ddiff(:).*A(:))-sum(wDT(:).*A(:)))/sum(Ddiff(:))));
invar2 = sprintf('\\textbf{invar} &= &%0.3f',1-var(Ddiff(:).*A(:)-wDT(:).*A(:))/var(Ddiff(:)));
text2 = strcat('\begin{tabular}{lcr}',invar2,'\\',bias2,'\end{tabular}');
text(256+80,128,text2,'interpreter','latex','rotation',90,'horizontalalignment','center')
axis off
caxis(cax)

subplot(3,3,9)
wDTA = waverec2(DT.*pass,sdiff,wavename);
imagesc(wDTA.*A)
axis image
hold on
plot(gx,gy,'k-')
% text and titles
title(sprintf('%d wavelets, masked',sum(DT.*pass~=0)))
bias3 = sprintf('\\textbf{bias} & = &%0.3e',abs((sum(Ddiff(:).*A(:))-sum(wDTA(:).*A(:)))/sum(Ddiff(:))));
invar3 = sprintf('\\textbf{invar} &= &%0.3f',1-var(Ddiff(:).*A(:)-wDTA(:).*A(:))/var(Ddiff(:)));
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
title(sprintf('buffered wavelet mask',sum(DT.*pass~=0)))
axis off
caxis(cax)

subplot(3,3,5)
[wA,s]=wavedec2(A,level,wavename);
AT = waverec2(wA.*pass,sdiff,wavename);
imagesc(cax(1)*AT)
axis image
hold on
plot(gx,gy,'w-')
% text and titles
title(sprintf('threshold by area',sum(DT.*pass~=0)))
bias4 = sprintf('\\textbf{bias} & = &%0.3e',abs((sum(A(:))-sum(AT(:)))/sum(A(:))));
invar4 = sprintf('\\textbf{invar} &= &%0.3f',1-var(A(:)-AT(:))/var(A(:)));
text4 = strcat('\begin{tabular}{lcr}',invar4,'\\',bias4,'\end{tabular}');
text(128,290,text4,'interpreter','latex','horizontalalignment','center')
axis off
caxis(cax)


colormap(bluewhitered(1000,1));


%%
% invariance curve for the classified image
level=8;
wavename = 'haar';
[wd,s]=wavedec2(A,level,wavename);

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
    wD{i} = waverec2(NC{i},sA,wavename);
    er(i) = immse(Ddiff,wD{i});
    r2(i) = 1-var(Ddiff(:)-wD{i}(:))/var(Ddiff(:));
    b(i) =  abs((sum(Ddiff(:))-sum(wD{i}(:)))/sum(Ddiff(:)));
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
