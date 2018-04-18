addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

load('Greenland60data');
load('ptsGL9')
load('im_tools9')
load('im_seqSH9.mat')

%%
in = 2;
for i = 1:size(D,3)
    X = D(:,:,i);
    [wd]=wavedec2(X,2,'haar');
    y(i) = wd(in);
end
figure(1)
clf
subplot(1,2,1)
wd([1:in-1,in+1:end]) = 0;
wD = waverec2(wd,s,'haar');
wd(in)
imagesc(wD(:,:));
imagesc(wD)
colormap(bluewhitered(1000,1));
axis image
hold on;
% Greenland
plot(gx,gy,'k-')
axis off

subplot(1,2,2)
plot(thedates,y,'o-','linewidth',1)
datetick
legend('wavelet coefficient value')
xlabel('year')
ylabel('kg per m^2')
%%
figure(2)
clf
imagesc(wD(:,:));
imagesc(wD)
colormap(bluewhitered(1000,1));
% colorbar
axis image
hold on;
% Greenland
plot(gx,gy,'k-')

%%
X = D(:,:,end);
[wd,s]=wavedec2(X,2,'haar');
% wd(abs(wd)<10*mean(abs(wd))) =0;

wd([1:in-1,in+1:end]) = 0;
wD = waverec2(wd,s,'haar');
wd(in)

figure(2)
clf
imagesc(X(:,:));
colormap(bluewhitered(1000,1));
colorbar
axis image
hold on;
% Greenland
plot(gx,gy,'k-')
axis off