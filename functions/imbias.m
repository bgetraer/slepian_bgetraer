function [ bias ] = imbias ( originalimage, testimage )
%IMINVAR Calculates the total energy loss (bias) between an original image
%and a test image. 
%   Image-to-image comparison of the total energy of the images. Return
%   valus can range from 0:+inf, where 0 would be no change of total pixel
%   value, and >1 would be a change in total pixel value greater than the
%   sum of the original image.
%
% INPUT
%   originalimage   the baseline image 
%   testimage       the image whose energy loss from the baseline is
%                       requested
% OUTPUT
%   bias           the image-to-image energy loss
%                       |(sum(O - T)/sum(O))|
%                  ex:
%                       bias(A,A*c) = 1-c
% last modified by: bgetraer@princeton.edu, 2/17/2019

% ensure images are comparable
if ~size(testimage) == size(originalimage)
    error('Images must be the same size!')
end

% bias = |energy(error)/energy(data)|
bias =  abs((sum(originalimage(:) - testimage(:)))/sum(originalimage(:)));
end