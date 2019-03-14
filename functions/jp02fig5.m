function jp02fig5(xprime,yprime,zprime,lond,latd,gx,gy,bx,by,Fx,Fy,n)
%JP02FIG5 Creates Figure 5 from "REGIONAL FORCING OF GREENLAND ICE LOSS 
%   2002-2017" Spring 2018 Junior Paper, Princeton Department of Geosciences
%
% last modified by: bgetraer@princeton.edu, 6/25/2018

defval('n',8);

% endpoints of the box
x1=Fx(lond(1,1),latd(1,1));
y1=Fy(lond(1,1),latd(1,1));
x2=Fx(lond(1,end),latd(1,end));
y2=Fy(lond(1,end),latd(1,end));
x3=Fx(lond(end,1),latd(end,1));
y3=Fy(lond(end,1),latd(end,1));
x4=Fx(lond(end,end),latd(end,end));
y4=Fy(lond(end,end),latd(end,end));

% Plot greenland and the outline of the box in Lat/Lon
figure(1)
clf 
hold on
greenland
% the outline
plot(lond(:,1),latd(:,1),'k-');
plot(lond(:,end),latd(:,end),'k-');
plot(lond(1,:),latd(1,:),'k-');
plot(lond(end,:),latd(end,:),'k-');
% the endpoints
plot(lond(1,1),latd(1,1),'o',...
    'MarkerF','r','MarkerE','k');
plot(lond(end,end),latd(end,end),'o',...
    'MarkerF','r','MarkerE','k');
plot(lond(1,end),latd(1,end),'o',...
    'MarkerF','r','MarkerE','k');
plot(lond(end,1),latd(end,1),'o',...
    'MarkerF','r','MarkerE','k');
title('Lat/Lon basis')

% Plot Greenland and the outline of the box in the 3d cartesian projection
figure(2)
clf 
subplot(1,2,1)
plotplmgreenland
hold on
plotcont([],[],4)
hold on
% the outline
plot3(xprime(:,1),yprime(:,1),zprime(:,1),'r-','linewidth',2);
plot3(xprime(:,end),yprime(:,end),zprime(:,end),'r-','linewidth',2);
plot3(xprime(1,:),yprime(1,:),zprime(1,:),'r-','linewidth',2);
plot3(xprime(end,:),yprime(end,:),zprime(end,:),'r-','linewidth',2);
% the endpoints
plot3(xprime(1,1),yprime(1,1),1,'o',...
    'MarkerF','r','MarkerE','k');
plot3(xprime(end,end),yprime(end,end),1,'o',...
    'MarkerF','r','MarkerE','k');
plot3(xprime(1,end),yprime(1,end),1,'o',...
    'MarkerF','r','MarkerE','k');
plot3(xprime(end,1),yprime(end,1),1,'o',...
    'MarkerF','r','MarkerE','k');
% formatting
axis([-0.17 0.55 -0.55 0.17 -1 1 ])
view(45,90)
colormap(bone)
title('Global basis')

% Plot Greenland and the outline of the box in the image projection
subplot(1,2,2)
axis square tight
set(gca,'ydir','reverse','box','on')
hold on

% endpoints
plot([x1,x2,x3,x4],[y1,y2,y3,y4],'o',...
    'MarkerF','r','MarkerE','k');
% Greenland
plot(gx,gy,'k-')
plot(bx,by,'k:')
linetik = [0 2.^(1:n)]+0.5;
tik = [0 2.^(5:n)]+0.5;
tiklab = [1 2.^(5:n)];
set(gca,'xtick',tik,'xticklabels',tiklab)
set(gca,'ytick',tik,'yticklabels',tiklab)
plotline1 = [linetik; linetik];
plotline2 = [ones(1,n+1); [linetik(2:end),linetik(end)]];
plot(plotline1,plotline2,'k')
plot(plotline2,plotline1,'k')
title('Image basis')

end

