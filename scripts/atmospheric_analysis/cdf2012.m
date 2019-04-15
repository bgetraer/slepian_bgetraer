%% SUMMER 2012 CDF
close all
run('compareTempMass.m')

%%

% mnth = [5 5];
% yr = [2003 2017];

startD = [6 2012];
endD = [9 2012];

mnth = [6 9];
yr = [2012 2012];

[meltJJA2012, massJJA2012] = meltdayvsmass(mnth,yr,massData.D(resamp,resamp,:),...
    massData.thedates,GRACEX,GRACEY,meltData.allmeltMap,meltData.alldates,pM.X,pM.Y);

figure(5)
clf
suptitle(sprintf('%d/%d-%d/%d',mnth(1),yr(1),mnth(2),yr(2)))
ax3 = subplot(2,4,1);
hold on
fill(pG.gx./resamprate+0.5,pG.gy./resamprate+0.5,[0.6 0.6 0.6],'EdgeColor','none')
set(gca,'ydir','reverse')
imagesc(massJJA2012.*indexICE,'AlphaData',(indexICE));
colormap(ax3, bluewhitered([],1))
plot(pG.bx./resamprate+0.5,pG.by./resamprate+0.5,'--k','linewidth',0.2)
axis image off
title('Mass Loss')

ax2 = subplot(2,4,2);
hold on
fill(pG.gx./resamprate+0.5,pG.gy./resamprate+0.5,[0.6 0.6 0.6],'EdgeColor','none')
set(gca,'ydir','reverse')
imagesc(meltJJA2012.*indexICE,'AlphaData',(indexICE));
colormap(ax2,parula)
plot(pG.bx./resamprate+0.5,pG.by./resamprate+0.5,'--k','linewidth',0.2)
axis image off
title('Melt Days')

ax2 = subplot(2,4,5:6);
hold on
plotGLbackground( pG,ice,[0.7,1,1] );
plotGLforeground(subREG,1:4,1);
set(gca,'ydir','reverse')
axis image off
subplot(1,2,2)

[~,x,y] = cdfMelt(meltJJA2012, massJJA2012, indexICE, subREG );


% yyaxis right
% 
% [a,i] = sort(meltJJA2012(:).*indexICE(:));
% b = distweight(:).*indexICEnan(:);
% b = b(i);
% 
% plot(a(~isnan(b)),smooth(b(~isnan(b)),20))
% axis tight square
legend('nw','ne','se','sw','location','southoutside')

%% EXAMINE SUMMER 2012
figure(6)
clf
split = [7,18];
for i=1:4
    Y = cumsum(y{i}(y{i}>0));
    X = x{i}(y{i}>0);
    
    subplot(1,2,1)
    hold on
    plot(X,Y,'linewidth',2);
    
    subplot(2,2,2)
    hold on
    if any(i==[1 3])
        plot(X(X<split(1)),Y(X<split(1)),'linewidth',2);
    else
        plot(X(X<split(2)),Y(X<split(2)),'linewidth',2);
    end
    axis tight
    
    subplot(2,2,4)
    hold on
    if any(i==[1 3])
        if i == 3
            plot(X,Y,'linewidth',2);
        else
            ysplit = Y(X>split(1));
            plot(X(X>split(1))-split(1),ysplit - ysplit(1),'linewidth',2);
        end
    else
        if i == 2
            plot(X,Y,'linewidth',2);
        else
            ysplit = Y(X>split(2));
            plot(X(X>split(2))-split(2),ysplit - ysplit(1),'linewidth',2);
        end
    end
    axis tight
end

subplot(1,2,1)
legend('nw','ne','se','sw','location','southoutside')
xlabel('number of melt days')
ylabel('kg per m^2 mass loss (cumulative)')

subplot(2,2,2)
xlabel('number of melt days')
ylabel('kg per m^2 mass loss (cumulative)')

subplot(2,2,4)
xlabel('number of melt days')
ylabel('kg per m^2 mass loss (cumulative)')
