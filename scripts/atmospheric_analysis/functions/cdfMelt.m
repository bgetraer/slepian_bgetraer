function [ hCDF,xsort,ysort ] = cdfMelt(meltdata, massdata, indexICE, subREG )
%CDFMELT Plots out a cumulative density function for comparing melt days to
% mass loss.

hold on
title('CDF')
hCDF = cell(1,4);
for i = 1:length(subREG.index)
    thisindex = indexICE.*subREG.index{i};
    size(thisindex(:))
    size(meltdata(:))
    x = meltdata(:).*thisindex(:);
    y = -massdata(:).*thisindex(:);
    [xsort{i},I] = sort(x);
    ysort{i} = y(I);
    Y = cumsum(ysort{i});
    hCDF = plot(xsort{i},Y,'-','linewidth',2);
end
xlabel('number of melt days')
ylabel('kg per m^2 mass loss (cumulative)')
end

