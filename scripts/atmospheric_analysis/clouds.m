
CLDTOT = load(fullfile(matDir,'CLDTOT2003-2017'),'data','t','spacelim');
CLDTOT.data = squeeze(CLDTOT.data);
%% CALCULATE AND PLOT MONTHLY ANOMALY
season = alldates;
season(month(alldates)>=12 | month(alldates)<=2) = 1;
season(month(alldates)>=3 & month(alldates)<=5) = 2;
season(month(alldates)>=6 & month(alldates)<=8) = 3;
season(month(alldates)>=9 & month(alldates)<=11) = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CALCULATE MONTHLY ANOMALY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ cloudAnom, avgCloud ] = anomMonth( CLDTOT.data, alldates);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT MONTHLY ANOMALY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
clf
cax = [0.35 0.85];

for m = 1:12
    subplot(3,4,m)
    imj = avgCloud(:,:,m).*GRISmerranan';
    
    thisimj = interp2(pM.X,pM.Y,imj',GRACEX,GRACEY);
    
    fill(pG.gx./resamprate+0.5,pG.gy./resamprate+0.5,[0.6 0.6 0.6],'EdgeColor','none')
    hold on
    plot(pG.bx./resamprate+0.5,pG.by./resamprate+0.5,':k','linewidth',0.2)
    set(gca,'ydir','reverse')
    
    imagesc(thisimj,'AlphaData',~isnan(GRISnan));
    
%     imPlot(imj,'jet',cax)
    hold on
%     plot(bx,by,'k--')
%     plot(gx,gy,'k')
    title(datestr(datenum(1,m,1),'mmmm'))
    caxis(cax)
%     colorbar('eastoutside')
    axis tight
    axis off
end

y = unique(year(alldates));

suptitle(sprintf('fraction of cloud cover, %i-%i',...
    min(y),max(y)))

axes('Position',[0.875 0.3 0.02 0.4])
axis off
caxis(cax)
cb = colorbar;
cb.Position = cb.Position + [0 0 0.01 0];
ylabel(cb,'fraction of total cloud cover','fontsize',10)

%% MAP OUT ALL YEARS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT MELT DAYS ANOMALY BY YEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
clf
cax =   [-0.15    0.15];

for i = 1:length(y)
    thisYear = y(i);
    thisYMelt = mean(cloudAnom(:,:,year(alldates)==thisYear & season==3),3);
    
    total(i) = sum(thisYMelt(:));
    
    subplot(3,5,i)
    
    imj = thisYMelt;
    
    thisimj = interp2(pM.X,pM.Y,imj',GRACEX,GRACEY);
    
    fill(pG.gx./resamprate+0.5,pG.gy./resamprate+0.5,[0.6 0.6 0.6],'EdgeColor','none')
    hold on
    plot(pG.bx./resamprate+0.5,pG.by./resamprate+0.5,':k','linewidth',0.2)
    set(gca,'ydir','reverse')
    
    imagesc(thisimj,'AlphaData',~isnan(GRISnan));
    
%     imPlot(imj,parula,cax)
    hold on
%     plot(bx,by,'k--')
%     plot(gx,gy,'k')
    title(num2str(thisYear))
    caxis(cax)
%     colorbar('eastoutside')
    axis tight  off
end

suptitle('Summer total cloud cover fraction anomaly by year, normalized by monthly mean')

axes('Position',[0.9 0.3 0.02 0.4])
axis off
caxis(cax)
cb = colorbar;
cb.Position = cb.Position + [0 0 0.01 0];
ylabel(cb,'fraction of cloud cover','fontsize',10)