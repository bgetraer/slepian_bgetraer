function [ lmcosi_diff ] = difflmcosi(lmcosi_data,c)
%DIFFLMCOSI subtract two lmcosi matrices of the same order
%   differences the coefficients, keeps the l's and m's, and scale by
%   constant c.
%
%INPUT
%   lmcosi_data     stack of lmcosi matrices (lmX6Xdates)
%   c               scaling constant (default: 1)
%
%OUTPUT
%   lmcosi_diff     the lmcosi matrix lmcosi1-lmcosi2
%
% Last modified by bgetraer@princeton.edu, 12/26/2017
defval('c',1)
lm = lmcosi_data(:,1:2,1);
for i = 1:size(lmcosi_data,3)-1
    lmcosi_data(:,1:2,i)
    % subtract and scale coefficients
    cosi_diff = c*(lmcosi_data(:,3:4,i)-lmcosi_data(:,3:4,i+1));
end
lmcosi_diff = [lm cosi_diff];
end

