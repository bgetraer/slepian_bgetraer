function [ invar ] = iminvar( originalimage, testimage )
%IMINVAR Calculates the image invariance between an original image and a
%test image. The original image is treated as the data, the test image as
%the model.
%   Pixel-to-pixel comparison returning an invariance value similar to a
%   coefficient of determination r^2. Values can run from -inf:1,
%   where 1 is perfect correlation, and values <0 reflect that the mean of
%   the original image would be a close predictor of image values than the
%   test image.
%
% INPUT
%   originalimage   the baseline image 
%   testimage       the image whose invariance from the baseline is
%                       requested
% OUTPUT
%   invar           the pixel-to-pixel invariance
%                       1 - var(T(:)-O(:))/var(T(:))
%                  ex:
%                       bias(A,A*(1-c)) = 1-c^2
% last modified by: bgetraer@princeton.edu, 2/17/2019

% ensure images are comparable
if ~size(testimage) == size(originalimage)
    error('Images must be the same size!')
end

% invariance = 1-var(error)/var(data)
invar = 1-var(originalimage(:)-testimage(:))/var(originalimage(:));

end

