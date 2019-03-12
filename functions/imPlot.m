function [  ] = imPlot( theimage, cmap,cax )
%IMDISPLAY Plot image rotating and using imagesc with independent colormap
%
% Last modified by bgetraer@princeton.edu 3/9/2019



defval('cmap','jet');
defval('cax',[prctile(theimage(:),3) prctile(theimage(:),98)]);

theimage = imrotate(theimage,90);
    
% plotting
imagesc(theimage,'AlphaData',~isnan(theimage));
hold on
% plot(bx,by,'w:')
% plot(contx,conty,'k-','linewidth',1)
    
% color
thisax = gca;
colormap(thisax, cmap)

caxis(cax)
% axis
% latlonaxis( xlim, ylim, thespacelim, 5 );

end

