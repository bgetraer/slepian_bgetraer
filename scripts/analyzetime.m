addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

order = 10;
load(strcat('ptsGL',num2str(order)))
load(strcat('im_tools',num2str(order)))
load(strcat('im_seqSH',num2str(order)))
load(strcat('thresh',num2str(order)))
load(strcat('pass',num2str(order),'lin2'))
%% image mask
% find the points inside of Greenland
a = inpolygon(xp(:),yp(:),bx,by);
A = 1*reshape(a,size(xp));
%% wavelet decomposition
level = 8;
wavename = 'haar';

wD = zeros(size(D,3),65536);

for i = 1:size(D,3)
    [wD(i,:), s] = wavedec2(D(:,:,i),level,wavename);
end

%% differencing
wDiff = wD - wD(1,:);
%% wavelet thresholding
index = threshpass.*pass;
wDT = wDiff.*index;

%% image reconstruction
DT = zeros(size(D));
for i = 1:size(D,3)
    DT(:,:,i) = waverec2(wDT(i,:),s,wavename);
end

%% mask the images
% the wavelet transformed images
DTmask = DT.*A;
% the original spherical harmonic images
Dmask = (D - D(:,:,1)).*A;

%% compare to originals
for i = 1:size(D,3)
    imr = Dmask(:,:,i);
    imt = DTmask(:,:,i);
    b(i) =  (sum(imr(:))-sum(imt(:)));%abs((sum(imr(:))-sum(imt(:))))/sum(imr(:)));
    invar(i) = 1-var(imr(:)-imt(:))/var(imr(:));
    totmasst(i) = sum(imt(:));
    totmassr(i) = sum(imr(:));
end
figure(10)
clf;hold on;
plot(thedates,invar)
ylim([0,1])

yyaxis right
plot(thedates,totmasst,'--k')
plot(thedates,totmassr,'-.k')
plot(thedates,b)
[m,f] = linear_m(thedates-thedates(1),b,1);
plot(thedates,f)
grid on
datetick
%% corrected by linear bias term

for i = 1:size(D,3)
    imr = Dmask(:,:,i) - f(i)./size(Dmask(:,:,i),1)^2;
    imt = DTmask(:,:,i);
    b(i) =  (sum(imr(:))-sum(imt(:)));%abs((sum(imr(:))-sum(imt(:))))/sum(imr(:)));
    invar(i) = 1-var(imr(:)-imt(:))/var(imr(:));
    totmassr(i) = sum(imr(:));
end
figure(11)
clf;hold on;
plot(thedates,invar)
ylim([0,1])

yyaxis right
plot(thedates,totmasst,'--k')
plot(thedates,totmassr,'-.k')
plot(thedates,b)

grid on
datetick

%% view image sequence
figure(1)
clf
subplot(1,2,1)
imagesc(Dmask(:,:,end));
axis image
cax = caxis;

subplot(1,2,2)
imagesc(DTmask(:,:,end));
axis image
caxis(cax);
colormap(bluewhitered(1000,1));
%% here is the "movie"

for i=1:size(DTmask,3)
    imagesc(DTmask(:,:,i));
%     caxis(cax)
    colormap(bluewhitered(1000,1));
    colorbar
    axis image
    hold on;axis off;
    % Greenland
    plot(gx,gy,'k-')
    title(sprintf('%s',datestr(thedates(i))))
    pause(0.1)
end

%% model each coeff
wDTM = zeros(size(wDT,1),size(wDT,2));
Mvar = zeros(sum(index),1);

cin = find(index==1);

for i=1:sum(index)
    [wDTM(:,cin(i)),Mvar(i)] = modCoeff(thedates,wDT(:,cin(i)));
end

ERR = 2*sum(sqrt(Mvar));

%% plot the model for a coeff if you want
figure(3)
clf 
hold on
errCOFF = cin(abnorm);
for j = errCOFF
    clf 
    hold on
    plot(thedates,wDT(:,j))
    plot(thedates,wDTM(:,j))
    datetick
    pause
end
% plot(cin,Mvar,'x')

%% model image reconstruction
DTM = zeros(size(D));

for i = 1:size(D,3)
    DTM(:,:,i)  = waverec2(wDTM(i,:),s,wavename);
end
%% mask the model images
% the wavelet transformed images
DTMmask = DTM.*A;

%% compare (with the corrected by linear bias term)

for i = 1:size(D,3)
    imr = Dmask(:,:,i) - f(i)./size(Dmask(:,:,i),1)^2;
    imt = DTMmask(:,:,i);
    b(i) =  (sum(imr(:))-sum(imt(:)));%abs((sum(imr(:))-sum(imt(:))))/sum(imr(:)));
    invar(i) = 1-var(imr(:)-imt(:))/var(imr(:));
    totmassr(i) = sum(imr(:));
    totmasst(i) = sum(imt(:));
end
figure(12)
clf;hold on;
plot(thedates,invar)
ylim([0,1])

yyaxis right
errorbar(thedates,totmassr,repmat(ERR,size(thedates)))
plot(thedates,totmasst,'-k','linewidth',1)
plot(thedates,totmassr,'-ok','linewidth',1,'markersize',2)

plot(thedates,b)
% plot(thedates,f)
grid on
datetick

figure(13)
clf
hold on
scatter(totmassr,totmasst)
plot(xlim,xlim)
axis([xlim xlim])
%% locate anomaly
% date of excursion
datein = year(thedates)>=2012 & year(thedates)<2014;

for i=cin
    y = wDT(:,i);
    sf = wDTM(:,i);
    abnorm(i==cin) = any(abs(y(datein)-sf(datein))>3*sqrt(Mvar(i==cin)));
    figure(20)
    clf
    hold on
    plot(thedates,y-sf)
    stdline(thedates,y-sf)
    datetick
    pause
end
    