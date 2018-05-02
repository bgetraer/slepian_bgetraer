addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

order = 10;
load(strcat('ptsGL',num2str(order)))
load(strcat('im_tools',num2str(order)))
load(strcat('im_endsSH',num2str(order)))
%%
level = 8;
% choose orthoganal wavelet for enery preservation
% we choose one that is blocky - we don't want to overfit our test case.
% choose fejer-korovkin for simplicity of representation (in minimum number
% of wavelets needed to cover the image, orthogonality,
% and improvement over haar in invariance.
wavename = 'haar';

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

level = 10;
wname = {'haar','fk4','bior6.8'};
ptl = [96:0.1:100];
N = 1:level;
    
figure(5)
clf
hold on
for j = 1:length(wname)
    Ddiff =  D(:,:,end)-D(:,:,1);
    [wdiff,sdiff]=wavedec2(Ddiff,level,wname{j});
    
    
    clear T er r2 wD b
    for i = 1:length(ptl)
        T(i) = prctile(wdiff,ptl(i));
        NC = wthcoef2('t',wdiff,sdiff,N,repmat(T(i),size(N)),'h');
        wD{i} = waverec2(NC,sdiff,wname{j});
        er(i) = immse(Ddiff,wD{i});
        r2(i) = 1-var(Ddiff-wD{i})/var(Ddiff);
    end
    loglog(ptl,r2,'linewidth',2)
end
%
legend(wname)

%%
figure(1)
imagesc(xp)
figure(2)
imagesc(yp)
hold on
plot(gx,gy,'k-')

a = inpolygon(xp(:),yp(:),gx,gy);
%%
A = reshape(a,size(xp));
figure(2)
clf
imagesc(A)
hold on
colormap(bone)
plot(gx,gy,'r','linewidth',2)