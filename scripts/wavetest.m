addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

load('ptsGL11')
load('im_tools11')
% load('im_seqSH9.mat')
load('im_endsSH11.mat')
%%
in = 3;
level = 10;
% choose orthoganal wavelet for enery preservation
% choose fejer-korovkin for 
wavename = 'fk4';
[phi,psi,xval] = wavefun(wavename);
figure(10)
subplot(2,1,1)
plot(xval,phi)
title('Scaling Function')
subplot(2,1,2)
plot(xval,psi)
title('Wavelet')

Ddiff =  D(:,:,end)-D(:,:,1);

figure(1)
clf
hold on
% imagesc(D(:,:,end));
% colormap(bluewhitered(1000,1));
% colorbar
axis image
set(gca,'ydir','reverse')

% Greenland
% plot(gx,gy,'k-')
[wdiff,sdiff]=wavedec2(Ddiff,level,wavename);
waveC = [];
l = zeros(level,1);
sz = zeros(level,2);
clear C 
for j = 1:level
    C{j} = appcoef2(wdiff,sdiff,wavename,j); % the wavelet coeff matrices
    sz(j,:) = size(C{j});    % save the sizes of each coeff matrix
    l(j) = sz(j,1)*sz(j,2); % save the lengths of each vector of coeffs
    waveC = [waveC C{j}(:)'];     % save the coeffs into a straight vector
    
    plot(gx,gy,'k-')
    imagesc(C{j})
    %     caxis([-0.2672    1.7162]*1E6)
    colormap(bluewhitered(1000,1));
end

figure(2)
clf
hold on
L = cumsum(l);
binwidth = 1;
for i = 1:length(l)
    if i==1
        x = linspace(0,binwidth,l(i));
        ordwaveC = sort(abs(waveC(1:l(i))));
    else
        x(end)
        x(end)+binwidth
        x = [x linspace(x(end),x(end)+binwidth,l(i))];
        ordwaveC = [ordwaveC sort(abs(waveC(L(i-1)+1:L(i))))];
    end
end

[sortwaveC, isort] = sort(abs(waveC));


plot(x,abs(waveC),'linewidth',0.5)
plot(x,ordwaveC,'linewidth',2)
plot(x,sortwaveC,'-.k','linewidth',2)
plot(xlim,[prctile(sortwaveC,99.5) prctile(sortwaveC,99.5)])

% set(gca,'xtick',a.^cumsum(l))
grid on
axis([0,level, ylim])

figure(3)
wD = waverec2(wdiff,sdiff,wavename);
imagesc(wD)
axis image
hold on
plot(gx,gy,'k-')
colormap(bluewhitered(1000,1));

% f = 1:9
% F = reshape(f,[3 3])
T = prctile(sortwaveC,90)
N = [1:level];
NC = wthcoef2('t',wdiff,sdiff,N,repmat(T,size(N)),'h');

figure(4)
wD = waverec2(NC,sdiff,wavename);
imagesc(wD)
axis image
hold on
plot(gx,gy,'k-')
colormap(bluewhitered(1000,1));


%%
ptl = [0:20:80 90:0.5:100];
clear T er r2 wD
for i = 1:length(ptl)
T(i) = prctile(sortwaveC,ptl(i));
N = [1:level];
NC = wthcoef2('t',wdiff,sdiff,N,repmat(T(i),size(N)),'h');
wD{i} = waverec2(NC,sdiff,wavename);

% er(i) = immse(Ddiff,wD);
r2(i) = 1-var(Ddiff-wD{i})/var(Ddiff);
end

%%
clf
hold on
yyaxis right
% plot(ptl,er,'linewidth',2)
yyaxis left
hold on
loglog(ptl,r2,'linewidth',2)

sum(sum(Ddiff))
sum(sum(wD{i}))
%%
index = 100:400;
level = 4;
wavename = 'fk4';
for in=index

% for i = 1:size(D,3)
%     X = D(:,:,i);
%     [wd,s]=wavedec2(X,level,wavename);
%     y(i) = wd(in);
% end
[wd,s]=wavedec2(Ddiff,level,wavename);

figure(1)
clf
subplot(1,2,1)
wd([1:in-1,in+1:end]) = 0;
wD = waverec2(wd,s,wavename);
% wd(in)
imagesc(wD(:,:));
imagesc(wD)
colormap(bluewhitered(1000,1));
axis image
hold on;
% Greenland
plot(gx,gy,'k-')
axis off
pause(0.01)
end
% subplot(1,2,2)
% plot(thedates,y,'o-','linewidth',1)
% datetick
% legend('wavelet coefficient value')
% xlabel('year')
% ylabel('kg per m^2')
% 
% 
% 
