function [ hBG ] = plotGLbackground( pG,ice,cmap,summit )
%PLOTGLBACKGROUND Plots a background of Greenland, buffer, and ice extent
% with an optional summit.

hBG = struct;
hBG.gl = fill(pG.gx,pG.gy,[0.6 0.6 0.6],'EdgeColor','none');
hBG.buff = plot(pG.bx,pG.by,'--k','linewidth',0.2);
hBG.ice = fill(ice.GRACE.X,ice.GRACE.Y,cmap(1,:));

if nargin > 3
    hBG.summit = plot3(summit(1),summit(2),2,'^r','markerfacecolor','r');
end

end

