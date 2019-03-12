addpath('/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/functions')
setworkspace('/Users/benjamingetraer/Documents/IndependentWork/SH_Workspace');
datadir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/datafiles';
addpath(datadir)

order = 10;
load(strcat('ptsGL'));%,num2str(order)))
load(strcat('im_tools'));%,num2str(order)))
load(strcat('im_seqSH'));%,num2str(order)))
load(strcat('thresh',num2str(order)))
load(strcat('pass',num2str(order),'lin2'))
%% image mask
% find the points inside of Greenland
a = inpolygon(xp(:),yp(:),bx,by);
A = 1*reshape(a,size(xp));
%% Mass transformation
% image basis area
imareaGL = polyarea(bx,by);

% image to global transform
im2glob = areaGL/imareaGL;

% kg to Gt
kg2Gt = 1E-12;

% both
imkg2globGt = im2glob*kg2Gt;
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
clf;
subplot(2,2,1:2)
hold on;
plot(thedates,invar,'linewidth',1)
ylim([0,1])
ylabel('invariance')

yyaxis right
plot(thedates,totmassr*imkg2globGt,'-k','linewidth',1)

plot(thedates,totmasst*imkg2globGt,'--k','linewidth',0.75)
plot(thedates,b*imkg2globGt,'linewidth',1)

[m,f] = linear_m(thedates-thedates(1),b,1);
plot(thedates,f*imkg2globGt,'linewidth',1)
grid on
legend('image invariance','total mass (full reconstruction)','total mass (thresholded reconstruction)','bias',sprintf('linear fit m=%0.0f Gt per year',m(2)),...
    'location','southwest')
title('Uncorrected for bias')
datetick
xlabel('year')
ylabel('mass (Gt)')
%% corrected by linear bias term

for i = 1:size(D,3)
    imr = Dmask(:,:,i);
    imt = DTmask(:,:,i)  + f(i)./size(Dmask(:,:,i),1)^2;
    b(i) =  (sum(imr(:))-sum(imt(:)));%abs((sum(imr(:))-sum(imt(:))))/sum(imr(:)));
    invar(i) = 1-var(imr(:)-imt(:))/var(imr(:));
    totmasst(i) = sum(imt(:));
end
subplot(2,2,3:4)
hold on;
plot(thedates,invar,'linewidth',1)
ylim([0,1])
ylabel('invariance')

yyaxis right
plot(thedates,totmassr*imkg2globGt,'-k','linewidth',1)

plot(thedates,totmasst*imkg2globGt,'--k','linewidth',0.75)
plot(thedates,b*imkg2globGt,'linewidth',1)

grid on
legend('image invariance','total mass (full reconstruction)','total mass (thresholded reconstruction)',sprintf('bias: \\sigma=%0.0f Gt',std(b*imkg2globGt)),'location','southwest')
title('Corrected for bias')
datetick
xlabel('year')
ylabel('mass (Gt)')

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
% errCOFF = cin(abnorm);
for j = cin(15)
    clf
    subplot(1,2,2)
    hold on
    plot(thedates,-wDT(:,j),'k','linewidth',0.75)
    plot(thedates,-wDTM(:,j),'-r','linewidth',1)
    datetick
    grid on
    lgd = legend('wavelet coefficient m_{\zeta=6,\gamma=5}',...
        'wavelet coefficient model');
    set(lgd,'fontsize',12,'location','southwest')
    subplot(1,2,1)
    testwD = zeros(1,65536);
    testwD(j) = -1;
    testD = waverec2(testwD,s,wavename);
    imagesc(testD)
    hold on
    plot(bx,by,'k:')
    plot(gx,gy,'k')
    colormap(bluewhitered(10,1))
    caxis([-1,1]*2E-2)
    
    text(96,13,'-','fontsize',30,'horizontalalignment','center')
    text(96,45,'+','fontsize',30,'horizontalalignment','center')
    
    set(gca,'xtick',[0.5 32.5 64.5 128.5 256.5],'xticklabels',[1 32 64 128 256])
    set(gca,'ytick',[0.5 32.5 64.5 128.5 256.5],'yticklabels',[1 32 64 128 256])
    plotline1 = [1.5 2.5 4.5 8.5 16.5 32.5 64.5 128.5 256.5; 1.5 2.5 4.5 8.5 16.5 32.5 64.5 128.5 256.5];
    plotline2 = [0 0 0 0 0 0 0 0 0; 2.5 4.5 8.5 16.5 32.5 64.5 128.5 256.5 256.5];
    plot(plotline1,plotline2,'k')
    plot(plotline2,plotline1,'k')
    % pause
end
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
    imr = Dmask(:,:,i);
    imt = DTMmask(:,:,i) + f(i)./size(Dmask(:,:,i),1)^2;
    b(i) =  (sum(imr(:))-sum(imt(:)));%abs((sum(imr(:))-sum(imt(:))))/sum(imr(:)));
    invar(i) = 1-var(imr(:)-imt(:))/var(imr(:));
    totmassr(i) = sum(imr(:));
    totmasst(i) = sum(imt(:));
end
figure(12)
clf;hold on;
% date of excursion
datex = minmax(thedates(year(thedates)>=2012 & year(thedates)<2014));
ylym = ylim;
h2012 = fill([flip(datex) datex],[ylym(1) ylym ylym(2)],[1 1 1]*0.9);

hinv = plot(thedates,invar,'b','linewidth',1);
ylim([0,1])
ylabel('invariance')

yyaxis right




he = errorbar(thedates,totmassr*imkg2globGt,repmat(ERR,size(thedates))*imkg2globGt,'color',[1 1 1]*0.4);
hr= plot(thedates,totmassr*imkg2globGt,'-k','linewidth',1,'markersize',2);

ht = plot(thedates,totmasst*imkg2globGt,'--r','linewidth',1);

hb = plot(thedates,b*imkg2globGt,'-','linewidth',1)
% plot(thedates,f)
grid on
datetick
legend([hinv he hr ht hb h2012],{
    'invariance','propogated variance from wavelet modeling',...
    'total mass (full reconstruction)',...
    'total mass (m_{\zeta\gamma} wavelet model)','residuals',...
    '2012--2014 anomaly'},...
    'location','southwest')
title('modeled wavelet reconstruction')% with 6 & 7 year periodic terms')



xlabel('year')
ylabel('mass (Gt)')



axes('position',[0.67 0.38 0.25 0.35])
hold on
scatter(totmassr,totmasst,50,'.')
r2 = corrcoef(totmassr,totmasst);
axis square tight
grid on
plot(xlim,xlim,'linewidth',1.5)
set(gca,'box','on','xticklabels','','yticklabels','')
xlabel('full reconstruction')
ylabel('m_{\zeta\gamma} wavelet model')
title(sprintf('cross plot r^2=%0.3f',r2(2).^2))

%% look for periodic structure within the residuals
% example:
degf=[];
psd=[];
for i = 1:1E4
    ndata = length(thedates);
    z = randn(ndata,1);
    window = hanning(ndata);
    [PSD,X] = periodogram(z,[],[],1/30);
    % figure(1)
    % clf
    % plot(X,PSD)
    % hold on
    % findpeaks(PSD,X);
    
    % degf = [degf var(z)];
    
    psd(:,i) = PSD;
end
%%
figure(1)
clf
plot(X,psd(:,1),'-')
%%
figure(2)
clf
subplot(1,3,1)
hold off
% what are the max values in each set?
psdmax = max(psd);

hst = histogram(psdmax,500);
hstval = hst.Values;
hstbe = hst.BinEdges;

xx = hstbe(1:end-1)+1/2*hst.BinWidth;
% plot(hstval)
cs = cumsum(hstval); % how many were already there
yy = 1-cs/cs(end);
plot([0 xx],[1 yy],'linewidth',1)

hold on
ptlpeak = prctile(psdmax,95);
plot([ptlpeak ptlpeak],ylim,'linewidth',1)

plot([ptlpeak ptlpeak],[interp1(xx,yy,ptlpeak) 0])
outside = yy(xx>ptlpeak);
fill([xx(xx>ptlpeak) ptlpeak ptlpeak ],[ outside 0 interp1(xx,yy,ptlpeak) ],'r')

xlabel('max peak height')
ylabel('set percentile')
lgd = legend('\begin{tabular}{l}probability of finding a \\ higher peak in a given set \end{tabular}',...
    sprintf('\\begin{tabular}{l} $5^{th}$ percentile \\\\ max peak height = %2.2f \\end{tabular}',ptlpeak));
set(lgd,'interpreter','latex')
view(-90, 90);
grid on
xlim([0 1600])
hold on
text(1500,0.5,...
    sprintf('\\bf{Bootstrap on %0.0d sets of \n normally distributed random samples}',1E4),...
    'horizontalalignment','center')
% mean(degf)

% compare to one of our residual power spectrum
for i=d(10)
    y = wDT(:,i);
    sf = wDTM(:,i);
    
    testy = (y-sf)/std(y-sf);
    
    [PSD,X,window] = periodogram(testy,[],[],1/30);
    
    subplot(1,3,2:3)
    plot(X,PSD,'-x')
    hold on
    plot(xlim,[ptlpeak ptlpeak],'linewidth',1)
    plot(freq,max(pk),'ok','markerfacecolor','k')

    %     findpeaks(PSD,X);
    xlabel('frequency (days)')
    ylabel('peak height')
    ylim([0 1600])
    title('Example of model residual power spectrum, m_{6,6}')
    legend('power spectrum of wavelet residuals',...
        sprintf('5^{th} percentile max peak height = %2.2f',ptlpeak),...
        sprintf('peak at period of %0.2f years',period),'location','southeast')
end

axes('position',[0.69 0.45 0.2 0.45])

[pk,loc]=findpeaks(PSD,X);

freq = loc(pk == max(pk));
yearlength = 365.2422; % in days
period = 1/freq/yearlength; %in years

[m,f,Rsq,power_m] = periodic_m(thedates,y-sf,1/freq);
plot(thedates,y-sf)
hold on
plot(thedates,f)

axis tight
datetick
legend('wavelet model residuals',sprintf('%0.2f year periodic model',period))


axes('position',[0.47 0.45 0.18 0.45])
 testwD = zeros(1,65536);
    testwD(22) = 1;
    testD = waverec2(testwD,s,wavename);
    imagesc(testD)
    axis image
    hold on
    plot(gx,gy,'k-')
    plot(bx,by,'k:')
    set(gca,'xtick',[0.5 32.5 64.5 128.5 256.5],'xticklabels',[1 32 64 128 256])
    set(gca,'ytick',[0.5 32.5 64.5 128.5 256.5],'yticklabels',[1 32 64 128 256])
    plotline1 = [8.5 16.5 32.5 64.5 128.5 256.5; 8.5 16.5 32.5 64.5 128.5 256.5];
    plotline2 = [ 1 1 1 1 1 1;  16.5 32.5 64.5 128.5 256.5 256.5];
    plot(plotline1,plotline2,'k')
    plot(plotline2,plotline1,'k')
    caxis([-1 1]*0.3E-1)
    colormap(bluewhitered(1000,1));
    axis square

%%
    P=[];
    A = [];
    M = [];
for i=cin
    y = wDT(:,i);
    sf = wDTM(:,i);
    
    testy = (y-sf)/std(y-sf);
    
    [PSD,X,window] = periodogram(testy,[],[],1/30);
    
%     figure(2)
%     clf
%     subplot(1,3,1:2)
%     plot(thedates,testy)
%     subplot(1,3,3)
%     plot(X,PSD)
%     hold on
%     plot(xlim,[ptlpeak ptlpeak])
%     
%     %     findpeaks(PSD,X);
%     xlabel('frequency')
%     ylabel('peak height')
%     ylim([0 1400])
%     pause
    
    if any(PSD>ptlpeak)
        [PSD,X,window] = periodogram(y-sf,[],[],1/30);
        
        [pk,loc]=findpeaks(PSD,X);
        
        freq = loc(pk == max(pk));
        yearlength = 365.2422; % in days
        
        period = 1/freq/yearlength; %in years
        % period = 6;
        % freq = 1/period/yearlength;
        
        [m,f,Rsq,power_m] = periodic_m(thedates,y-sf,1/freq);
        
        
        P(i) = period;
        A(i) = max(pk);
        Rsq;
        
        M(:,i) = m;
        %     figure(21)
        %     clf
        %     plot(X,PSD)
        %     hold on
        %     findpeaks(PSD,X);
        
        %     figure(20)
        %     clf
        %     hold on
        %     title(strcat('wavelet',num2str(i),' : period=',num2str(period),' : r^2=',num2str(Rsq)))
        %     plot(thedates,y-sf)
        %     stdline(thedates,y-sf);
        %     plot(thedates,f,'--r');
        %     ylabel('model residuals')
        %     datetick
        %     pause
    end
end


%% check out the period structure
figure(23)
    clf
uP = unique(P);
cax = [-1 1]*10E7;

for k=1:length(uP)
    whereP = find(P==uP(end+1-k));
    
    testwD = zeros(1,65536);
    testwD(whereP) = A(whereP);
    
    testD = waverec2(testwD,s,wavename);
    
    
    if k~=6,    
        subfig = [[k:1+k]+k-1 [k:1+k]+k-1+6];
        if k>3, subfig = [[k:1+k]+k-1+18 [k:1+k]+k-1+24];
        end
        subplot(8,6,subfig)
        imagesc(testD)
        axis image
        hold on
        plot(gx,gy,'k-')
        plot(bx,by,'k:')
        set(gca,'xtick',[0.5 32.5 64.5 128.5 256.5],'xticklabels',[1 32 64 128 256])
        set(gca,'ytick',[0.5 32.5 64.5 128.5 256.5],'yticklabels',[1 32 64 128 256])
        plotline1 = [8.5 16.5 32.5 64.5 128.5 256.5; 8.5 16.5 32.5 64.5 128.5 256.5];
        plotline2 = [ 1 1 1 1 1 1;  16.5 32.5 64.5 128.5 256.5 256.5];
        plot(plotline1,plotline2,'k')
        plot(plotline2,plotline1,'k')
        ptext = sprintf('Period = %0.2f years',uP(end+1-k));
        ntext = sprintf('%0.0f Wavelet(s)',length(whereP));
        title(strcat('\begin{tabular}{c}',ntext,' \\ ',ptext,' \end{tabular}'),...
            'interpreter','latex')
        
        caxis(cax);
        
        subfig = [[k:1+k]+k-1+12]% [k:1+k]+k-1+18]
        if k>3, subfig = [[k:1+k]+k-1+30]% [k:1+k]+k-1+36]
        end
        
        whereP = find(P==uP(end+1-k));
        phase{k} = atan(M(1,whereP)./M(2,whereP));
        
        subplot(8,6,subfig)

        hold on
        xxx = 0:1:yearlength*15;
        yyy = sin(phase{k} + 2*pi/uP(end+1-k)/yearlength.*thedates');
        plot(thedates,yyy,'k')
        datetick
       axis tight 
       yticklabels '';grid on
        
        
        
        whereP = find(P==uP(end+1-k));
%         phase{k} = atan(M(1,whereP)./M(2,whereP));
        
        subplot(8,6,[41:42])
        title('All periods and phases')
        hold on
        c = flip(hsv(5));
        yyy = sin(phase{k} + 2*pi/uP(end+1-k)/yearlength.*thedates');
        plot(thedates,yyy,'color',c(k,:))
       datetick
       axis tight 
       yticklabels '';grid on
    end
    
    if k==6,   
        subplot(8,6,29:30)
        axis off
        caxis(cax);
        cb = colorbar;
        ylabel(cb,'\begin{tabular}{l} amplitude of \\ periodic structure \end{tabular}',...
            'interpreter','latex','fontsize',12,'rotation',0)
        set(cb,'location','north')%'position', [0.67 0.2 0.0183 0.3417],'AxisLocation', 'in')
    end
    
        colormap(bluewhitered([],1));
end  
%% 
for k=1%:length(uP)-1
    
    whereP = find(P==uP(end+1-k));
    phase{k} = atan(M(1,whereP)./M(2,whereP));
    
    figure(4)
    clf 
    hold on
    xxx = 0:0.1:4*pi;
    yyy = sin(phase{k} + xxx')
    plot(xxx,yyy)
%     
%     for i = 1:length(whereP)
%         figure(3)
%         clf
%         y = wDT(:,whereP(i));
%         sf = wDTM(:,whereP(i));
%         plot(thedates,y-sf)
%         hold on
%         amp = sqrt(M(1,whereP(i))^2 + M(2,whereP(i))^2 );
%         plot(thedates, amp*sin(mean(phase{k}) + 2*pi/uP(end+1-k)/yearlength * thedates),'r')
%         pause
%     end
end
%% re-model each coeff with the altered model
wDTM = zeros(size(wDT,1),size(wDT,2));
Mvar = zeros(sum(index),1);

cin = find(index==1);

for i=1:sum(index)
    [wDTM(:,cin(i)),Mvar(i)] = modCoeff(thedates,wDT(:,cin(i)),1);
end

ERR = 2*sum(sqrt(Mvar));

%% model image reconstruction
DTM = zeros(size(D));

for i = 1:size(D,3)
    DTM(:,:,i)  = waverec2(wDTM(i,:),s,wavename);
end
%% mask the model images
% find the points inside of Greenland
a = inpolygon(xp(:),yp(:),bx,by);
A = 1*reshape(a,size(xp));
% the wavelet transformed images
DTMmask = DTM.*A;

%% compare (with the corrected by linear bias term)

for i = 1:size(D,3)
    imr = Dmask(:,:,i);
    imt = DTMmask(:,:,i) + f(i)./size(Dmask(:,:,i),1)^2;
    b(i) =  (sum(imr(:))-sum(imt(:)));%abs((sum(imr(:))-sum(imt(:))))/sum(imr(:)));
    invar(i) = 1-var(imr(:)-imt(:))/var(imr(:));
    totmassr(i) = sum(imr(:));
    totmasst(i) = sum(imt(:));
end

figure(12)
clf;hold on;
% date of excursion
datex = minmax(thedates(year(thedates)>=2012 & year(thedates)<2014));
ylym = ylim;
h2012 = fill([flip(datex) datex],[ylym(1) ylym ylym(2)],[1 1 1]*0.9);

hinv = plot(thedates,invar,'b','linewidth',1);
ylim([0,1])
ylabel('invariance')

yyaxis right




he = errorbar(thedates,totmassr*imkg2globGt,repmat(ERR,size(thedates))*imkg2globGt,'color',[1 1 1]*0.4);
hr= plot(thedates,totmassr*imkg2globGt,'-k','linewidth',1,'markersize',2);

ht = plot(thedates,totmasst*imkg2globGt,'--r','linewidth',1);

hb = plot(thedates,b*imkg2globGt,'-','linewidth',1)
% plot(thedates,f)
grid on
datetick
legend([hinv he hr ht hb h2012],{
    'invariance','propogated variance from wavelet modeling',...
    'total mass (full reconstruction)',...
    'total mass (m_{\zeta\gamma} wavelet model)','residuals',...
    '2012--2014 anomaly'},...
    'location','southwest')
title('modeled wavelet reconstruction')% with 6 & 7 year periodic terms')



xlabel('year')
ylabel('mass (Gt)')



axes('position',[0.67 0.38 0.25 0.35])
hold on
scatter(totmassr,totmasst,50,'.')
r2 = corrcoef(totmassr,totmasst);
axis square tight
grid on
plot(xlim,xlim,'linewidth',1.5)
set(gca,'box','on','xticklabels','','yticklabels','')
xlabel('full reconstruction')
ylabel('m_{\zeta\gamma} wavelet model')
title(sprintf('cross plot r^2=%0.3f',r2(2).^2))



%



% 
% figure(12)
% clf;hold on;
% plot(thedates,invar)
% ylim([0,1])
% 
% yyaxis right
% errorbar(thedates,totmassr*imkg2globGt,repmat(ERR,size(thedates))*imkg2globGt)
% plot(thedates,totmasst*imkg2globGt,'--k','linewidth',0.5)
% plot(thedates,totmassr*imkg2globGt,'-ok','linewidth',1,'markersize',2)
% 
% plot(thedates,b*imkg2globGt)
% % plot(thedates,f)
% grid on
% datetick
% legend('invariance','propogated variance from wavelet modeling','total mass from wavelet modeling','total mass SH','difference',...
%     'location','southwest')
% title('modeled wavelet reconstruction with 6 & 7 year periodic terms')
% 
% figure(13)
% clf
% hold on
% scatter(totmassr,totmasst)
% r2 = corrcoef(totmassr,totmasst);
% plot(xlim,xlim)
% axis([xlim xlim])
% xlabel('total mass full reconstruction')
% ylabel('total mass modeled wavelets')
% title(strcat('cross plot model2 r^2=',num2str(r2(2).^2)))

%% save data
filename = sprintf('wavetimeseries%d',order);
save(fullfile(datadir,filename),'thedates','totmasst','totmassr',...
    'imkg2globGt','ERR','b','wDT','wDTM','s','wavename','cin','index','A','f')
