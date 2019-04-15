function [ hFG ] = plotGLforeground(subREG, region,label,color )
%PLOTGLFOREGROUND Plots a foreground of Greenland subregions
defval('label',0)
defval('color','k')

hFG = cell(1,length(region));

regionlabel = {'NW','NE','SE','SW'};
regioncenterX = [105 140 155 120 ];
regioncenterY = [90 80 130 145];

for reg = region
       hFG{reg} = fill3(subREG.GRACE.XBUF{reg},subREG.GRACE.YBUF{reg},...
            repmat(2,size(subREG.GRACE.XBUF{reg})),'w','FaceAlpha',0.5,...
            'edgecolor',color{reg},'linewidth',1);
        if label
            text(regioncenterX(reg),regioncenterY(reg),...
                3,regionlabel{reg});
        end
end



end
