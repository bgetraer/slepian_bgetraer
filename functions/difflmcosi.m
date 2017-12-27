function [ lmcosi_diff ] = difflmcosi(lmcosi1,lmcosi2)
%difflmcosi subtract two lmcosi matrices of the same order
%   differences the coefficients, keeps the l's and m's, and keeps the
%   order of the smaller matrix
%
%INPUT
%   lmcosi1     an lmcosi matrix
%   lmcosi2     another lmcosi matrix
%OUTPUT
%   lmcosi_diff     the lmcosi matrix lmcosi1-lmcosi2
%
% Last modified by bgetraer@princeton.edu, 12/26/2017
    
dimension = min([size(lmcosi1,1) size(lmcosi2,1)]);
lmcosi_diff(:,1:2) = lmcosi1(1:dimension,1:2);
lmcosi_diff(:,3:4) = lmcosi1(1:dimension,3:4)-lmcosi2(1:dimension,3:4);
end

