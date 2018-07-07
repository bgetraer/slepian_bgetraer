%% Locate the 2012 to 2014 anomaly in space
addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

%%

load(fullfile(datadir,sprintf('wavetimeseries%d',order)));

% date of excursion
datein = year(thedates)>=2012 & year(thedates)<2014;

% only check out the dates we want
x = thedates(datein);
cM = wDTM(datein,:);
cW = wDT(datein,:);
F = f(datein);

t = totmasst(datein);
r = totmassr(datein);

r2_ref = corrcoef(t,r);
R2_ref = r2_ref(2)^2;

for i=cin
    % the modeled coefficients
    cm  = cM;
    % replace with one of the non-modeled coefficients
    cm(:,i) = cW(:,i);
    
    % model image reconstruction
    cmim = zeros([256 256 sum(datein)]);
    
    for j = 1:size(cm,1)
        cmim(:,:,j)  = waverec2(cm(j,:),s,wavename);
    end
    % mask the model images
    cmimmask = cmim.*A;
    
    clear totmasst 
    for k = 1:size(cm,1)
        imt = cmimmask(:,:,k) + F(k)./size(A,1)^2;
%         b(i) =  (sum(imr(:))-sum(imt(:)));%abs((sum(imr(:))-sum(imt(:))))/sum(imr(:)));
%         invar(i) = 1-var(imr(:)-imt(:))/var(imr(:));
        totmasst(k) = sum(imt(:));
    end
    
    
    r = corrcoef(totmassr(datein),totmasst);
    R2(i) = r(2)^2;
    
%     abnorm(i==cin) = any(abs(y(datein)-sf(datein))>3*sqrt(Mvar(i==cin)));
%     figure(20)
%     clf  
%     hold on
%     title(num2str(i))
%     plot(thedates,y-sf)
%     stdline(thedates,y-sf);
%     datetick
%     pause
end
%%
figure(1)
clf;hold on;
errorbar(x,totmassr(datein)*imkg2globGt,repmat(ERR,size(x))*imkg2globGt)
plot(x,totmassr(datein)*imkg2globGt,'-ok','linewidth',1,'markersize',2)
plot(x,totmasst*imkg2globGt,'--k','linewidth',0.5)

grid on
datetick
% legend('invariance','propogated variance from wavelet modeling','total mass from wavelet modeling','total mass SH','difference',...
%     'location','southwest')
title('Reconstruction of the 2012-2014 Deviation')


% the improvements
impr = R2(cin)-R2_ref;

% figure(2)
% plot([0 impr])
% MPP = prctile(abs(impr),95);
% [~,loc] = findpeaks([0 impr],'MinPeakProminence',MPP);
% 
% loc = find(impr>0.022);
% 
% testwD = zeros(1,65536);
% testwD(cin(loc)) = cW(8,cin(loc));
% testD = waverec2(testwD,s,wavename);

% the modeled coefficients
    cm  = cM;
    % replace the desired non-modeled coefficients
    cm(:,cin(loc)) = cW(:,cin(loc));
    
    % model image reconstruction
    cmim = zeros([256 256 sum(datein)]);
    
    for j = 1:size(cm,1)
        cmim(:,:,j)  = waverec2(cm(j,:),s,wavename);
    end
    % mask the model images
    cmimmask = cmim.*A;
    
    for k = 1:size(cm,1)
        imt = cmimmask(:,:,k) + F(k)./size(A,1)^2;
%         b(i) =  (sum(imr(:))-sum(imt(:)));%abs((sum(imr(:))-sum(imt(:))))/sum(imr(:)));
%         invar(i) = 1-var(imr(:)-imt(:))/var(imr(:));
        totmasst2(k) = sum(imt(:));
    end
    
    
    r = corrcoef(totmassr(datein),totmasst2);
    r2 = r(2)^2;
    r = corrcoef(totmassr(datein),totmasst);
    r1 = r(2)^2;
    
 
    
figure(1)
hold on
plot(x,totmasst2*imkg2globGt,'-r','linewidth',1)
datetick('x', 'mmm yy')
xlabel('years')
ylabel('Gt since January 2003')
legend('propogated variance from wavelet modeling',...
    'total mass (full reconstruction)',sprintf('total mass (m_{\\zeta\\gamma} altered wavelet model): R^2=%0.2f',r1),...
    sprintf('total mass (altered m_{\\zeta\\gamma} model with 4 complete m_{\\zeta\\gamma} returned): R^2=%0.2f',r2),...
    'location','southoutside')
%%
figure(1)

axes('position',[0.6 0.6 0.25 0.35])
hold on

imagesc(testD)
set(gca,'ydir','reverse')
axis image 
box on

% hold on
hold on
        plot(gx,gy,'k-')
        plot(bx,by,'k:')
        set(gca,'xtick',[0.5 32.5 64.5 128.5 256.5],'xticklabels',[1 32 64 128 256])
        set(gca,'ytick',[0.5 32.5 64.5 128.5 256.5],'yticklabels',[1 32 64 128 256])
        plotline1 = [8.5 16.5 32.5 64.5 128.5 256.5; 8.5 16.5 32.5 64.5 128.5 256.5];
        plotline2 = [ 1 1 1 1 1 1;  16.5 32.5 64.5 128.5 256.5 256.5];
        plot(plotline1,plotline2,'k')
        plot(plotline2,plotline1,'k')
% pause
caxis([-1000 500])
colormap(bluewhitered([],1))
cb = colorbar;
cb.Position =  [0.86   0.63    0.0218    0.3];
ylabel(cb,'kg per m^2')
% figure(101)
% plot(x,cW(:,cin(loc))')

%%


