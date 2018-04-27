addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

order = 9;
load(strcat('ptsGL',num2str(order)))
load(strcat('im_tools',num2str(order)))
load(strcat('im_endsSH',num2str(order)))

close all
%%
% get the greenland buffer
buff=greenland(10,0.5);
bx = Fx(buff(:,1),buff(:,2));
by = Fy(buff(:,1),buff(:,2));

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

%%
% DETERMINE WHICH WAVELETS PASS AN AREA THRESHOLD WITHIN A BUFFERED
% GREENLAND
[wd,s]=wavedec2(A,level,wavename);
wd = wd*0+1;
pass = zeros(size(wd));
areathresh = 1/3;

for in=44
    fprintf('loop %d \n',in)
    % for i = 1:size(D,3)
    %     X = D(:,:,i);
    %     [wd,s]=wavedec2(X,level,wavename);
    %     y(i) = wd(in);
    % end

    figure(6)
    clf
    thiswd = wd;
    thiswd([1:in-1,in+1:end]) = 0;
    wD = waverec2(thiswd,s,wavename);
    xw = xp(wD(:)~=0);
    yw = yp(wD(:)~=0);
    pw = inpolygon(xw(:),yw(:),bx,by);
    aw = sum(wD(:)~=0);
    ainp = sum(pw(:)~=0);
    pass(in) = (1-(aw-ainp)/aw > areathresh);

    imagesc(~wD)
    colormap(bone);
    axis image
    hold on;
    % Greenland
    plot(gx,gy,'k-')
    axis off
    pause(0.1)
end

% sum(wD(:)~=0)
% subplot(1,2,2)
% plot(thedates,y,'o-','linewidth',1)
% datetick
% legend('wavelet coefficient value')
% xlabel('year')
% ylabel('kg per m^2')

% level 8
in8 = 1:4*l(8);
pass(in8)=1;
% level 7
in7 = in8(end)+ (1:3*l(7));
pass(in7)=1;
% level 6
in6 = in7(end)+ (1:3*l(6));
pass(in6)=0;
% level 6
in6 = in7(end)+ (1:3*l(6));
pass(in6)=0;
% level 6
in6 = in7(end)+ (1:3*l(6));
pass(in6)=0;


filename = sprintf('pass%d%d',order,areathresh);
save(fullfile(datadir,filename),'pass')
%%
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



%%
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
ptile = ptl(qc);
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
