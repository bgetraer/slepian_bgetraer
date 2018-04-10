addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

load('Greenland60data');
load('ptsGL9')
load('im_tools9')
load('im_seq9.mat')

%%
X = D(:,:,end);

[wd,s]=wavedec2(X,2,'haar');
wd(abs(wd)<mean(abs(wd))) =0;
wD = waverec2(wd,s,'haar');

figure(2)
clf
imagesc(wD(:,:));
    imagesc(wD)
    colormap(bluewhitered(1000,0));
    colorbar
    axis image
    hold on;
    % Greenland
    plot(gx,gy,'k-')
