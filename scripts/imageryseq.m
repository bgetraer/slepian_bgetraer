addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

load('Greenland60data');
load('ptsGL9')
load('im_tools9')
%% find peaks in Greenland signal
figure(1)
clf
plot(thedates,ESTtotal)
findpeaks(ESTtotal)
% store index
[pk,loc]=findpeaks(ESTtotal);

%% Get the LMCOSI coefficients for every image
shcoffs = G(:,1:20)*slepcoffs(:,1:20)';

% Create blank LMCOSI matrix
L=60;
[~,~,~,blank,~,~,~,~,~,ronm]=addmon(L);
lmcosi_mat = zeros([size(blank) 1]) + blank;

% Create the coefficient blanks
cosi=blank(:,3:4);

% creat the full signal matrix
full_lmcosi = zeros([size(lmcosi_mat), size(shcoffs,2)]);

for i=1:size(shcoffs,2)
    % clear dummy matrices
    this_mat = lmcosi_mat;
    this_blank = cosi;
    % grab the coefficients of an alpha eigentaper and
    % re-index them in lmcosi format
    this_blank(ronm)=shcoffs(:,i);
    % Add them to the full matrix
    this_mat(:,3:4,1)=this_blank;
    
    % add to the full signal matrix
    full_lmcosi(:,:,i) = this_mat;
end

%% Evaluate the signals on the grid
% blank for the images
shp = size(lond);
D = zeros([shp,size(shcoffs,2)]);
Dim = zeros(size(D));
% this is the time consuming bit... plm2xyz takes a while
for i=1:10
    fprintf('i=%d of %d \nnow starting plm2xyz\n',i,size(shcoffs,2))
    [data]=plm2xyz(full_lmcosi(:,:,i),latd(:),lond(:)); % vector of solutions
    D(:,:,i) = reshape(data,shp);  % put the vector back into matrix form
end
%% Save images
filename = 'im_seq9';
save(fullfile(datadir,filename),'D','full_lmcosi')
%%

% Greenland lat lon
gxy = greenland(10);
% Greenland xy in the image basis
gx = Fx(gxy(:,1),gxy(:,2));          
gy = Fy(gxy(:,1),gxy(:,2));

figure(2)
for i=1:size(D,3)
    clf
    imagesc(D(:,:,i));
    colormap(bluewhitered(1000,0));
    colorbar
    axis image
    hold on;axis off;
    % Greenland
    plot(gx,gy,'k-')
    pause(0.1)
end
%%
month(thedates)

% how mass has changed
dp = diff(pk);
mean(dp)
std(dp)
