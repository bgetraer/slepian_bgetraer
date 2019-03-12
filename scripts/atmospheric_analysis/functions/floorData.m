function [ datamatrix2 ] = floorData( datamatrix3 )
%FLOORDATA Takes an lxmxn matrix and returns an lxm matrix with the first
% real values available in the n layers.
%   Replaces NaN values in the lowest layer of the stack with real values
%   from the next layer up, iteratively until all NaN values are elimated
%   or there are no more layers.
%
% KNOWN BUGS
%   Data at the edges of layer transitions is not continuous
%
% INPUT
%   datamatrix3     lxmxn matrix of data
% OUTPUT
%   datamatrix2     lxm matrix of data
%
% Last modified by bgetraer@princeton.edu 3/1/2019

lvl = 1;
datamatrix2 = datamatrix3(:,:,lvl);

while any(isnan(datamatrix2(:))) && lvl <= size(datamatrix3,3)    
    slicenext = datamatrix3(:,:,lvl+1);
    msk = isnan(datamatrix2);
    datamatrix2(msk) = slicenext(msk);
    lvl=lvl+1;
end
end

