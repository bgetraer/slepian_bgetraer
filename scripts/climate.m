%% HADLEY SST DATA from 
%       https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html

% load files
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles/HadleyCentre';
addpath(datadir)

yr = 2004:2017;

for theyr = yr
    filename = sprintf('HadISST1_SST_%04d.txt',2004);

sst_had = importdata(fullfile(datadir,filename));
sstdata = sst_had.data(1:end-1,:);
%% SST anomaly
latsst = linspace(0.5,179.5,180);
lonsst = [linspace(0.5,359.5,360)];
figure(1)
clf
plot(lonsst)
%%
[lonmat,latmat] = ndgrid(lonsst,latsst);
latmat = latmat';
lonmat = lonmat';
%%
figure(1)
clf
gxy = greenland([],0.5);

% imagesc(latmat)
% scatter(lonmat(:),latmat(:),[],sstdata(:))
axis image
hold on
plot(gxy(:,1)-180,90-gxy(:,2),'r','linewidth',1)
hold on
%%

plot(lond-180,90-latd,'.')
%%
vq = griddata(lonmat,latmat,sstdata,lond-180,90-latd);

figure(1)
clf
imagesc(vq)
hold on
plot(gx,gy,'r','linewidth',1)
plot(bx,by,'r','linewidth',2)