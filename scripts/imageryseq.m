addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

load('Greenland60data');
load('ptsGL10')
load('im_tools10')
%% find peaks in Greenland signal
figure(1)
clf
plot(thedates,ESTtotal)
findpeaks(ESTtotal)
% store index
[pk,loc]=findpeaks(ESTtotal);

%% Get the LMCOSI coefficients for every image
% shcoffs = G(:,1:20)*slepcoffs(:,1:20)';
% 
% % Create blank LMCOSI matrix
L=60;
% [~,~,~,blank,~,~,~,~,~,ronm]=addmon(L);
% lmcosi_mat = zeros([size(blank) 1]) + blank;
% 
% % Create the coefficient blanks
% cosi=blank(:,3:4);
% 
% % creat the full signal matrix
% full_lmcosi = zeros([size(lmcosi_mat), size(shcoffs,2)]);
% 
% for i=1:size(shcoffs,2)
%     % clear dummy matrices
%     this_mat = lmcosi_mat;
%     this_blank = cosi;
%     % grab the coefficients of an alpha eigentaper and
%     % re-index them in lmcosi format
%     this_blank(ronm)=shcoffs(:,i);
%     % Add them to the full matrix
%     this_mat(:,3:4,1)=this_blank;
%     
%     % add to the full signal matrix
%     full_lmcosi(:,:,i) = this_mat;
% end

[potcoffs,cal_errors,thedates]=grace2plmt_inprogress('CSR','RL05','SD',0,L);

validrange = monthnum(1,2003,thedates):length(thedates); % crop off the first 7 months of untrustworthy data: 
potcoffs = potcoffs(validrange,:,:);
thedates = thedates(validrange);
%% Evaluate the signals on the grid
clear alphavar alphavarall blank cosi ESTresid ESTsignal ESTtotal ESTtotalresid ...
    ftests G lmcosi_mat lmcosi_sum slepcoffs signal 
% blank for the images
shp = size(lond);
filename = 'im_endsSH10';
% load(fullfile(datadir,filename))
D = zeros([shp,2]);
% this is the time consuming bit... plm2xyz takes a while
for i=[1 size(thedates,2)]
    fprintf('i=%d of %d \nnow starting plm2xyz\n',i,size(thedates,2))
    [data]=plm2xyz(squeeze(potcoffs(i,:,:)),latd(:),lond(:)); % vector of solutions
    D(:,:,i) = reshape(data,shp);  % put the vector back into matrix form
    % Save images
    save(fullfile(datadir,filename),'D','thedates')
end
%%

% Greenland lat lon
gxy = greenland(10);
% Greenland xy in the image basis
gx = Fx(gxy(:,1),gxy(:,2));          
gy = Fy(gxy(:,1),gxy(:,2));

Dprime = D-D(:,:,1);

figure(2)
for i=1:size(D,3)
    clf
    imagesc(Dprime(:,:,i));
    caxis([-3.8344    1.5619]*1E3)
    colormap(bluewhitered(1000,1));
    colorbar
    axis image
    hold on;axis off;
    % Greenland
    plot(gx,gy,'k-')
    title(sprintf('%s',datestr(thedates(i))))
    pause(0.1)
end
%%
month(thedates)

% how mass has changed
dp = diff(pk);
mean(dp)
std(dp)
