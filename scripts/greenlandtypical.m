load(fullfile(datadir,'Greenland60data'))

figure(1)

clf
allyears = unique(year(thedates));

for i = 1:14
    thisyear = allyears(i);
    seasonaldiff(i,:) = ESTsignal(monthnum(1,thisyear+1,thedates),1:20) - ESTsignal(monthnum(1,thisyear,thedates),1:20);
end

meandiff = sum(seasonaldiff,1)/14;



grabanalpha(G(:,1:20)*meandiff',[],60,1,2)
cmap = 'bluewhitered(20,1)';
colormap(eval(cmap))
view(45,90)
map_r = 0.4;
axis([-0 map_r -map_r 0])
axis off

buff=greenland(10,0.5);
gl=greenland(10,0);
XY1 = [buff;gl;];
XY2 = gl;
[x1,y1,z1] = sph2cart(XY1(:,1)*pi/180,XY1(:,2)*pi/180,1); 
[x2,y2,z2] = sph2cart(XY2(:,1)*pi/180,XY2(:,2)*pi/180,1); 
hold on
plot3(x2,y2,z2,'-k','linewidth',2); plot3(x1,y1,z1,'-k','linewidth',2);
patch(x1,y1,z1,[1,1,1]);
view(45,90)
axis off