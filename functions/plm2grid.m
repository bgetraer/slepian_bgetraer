function [ D ] = plm2grid( coffs, thedates, lat, lon )
%PLM2GRID Creates a 3d matrix of images (x,y,time) by evaluating the given
%   coefficient timeseriec upon a given lat lon grid. This is fairly time
%   and memory consuming, so the datafile is saved after every date is
%   done, and progress is regularly printed.
%
% INPUT:
%   coffs       the PLM coefficients from GRACE2PLMT or wherever
%   thedates    the corresponding dates 
%   lat         the latitude matrix of your grid
%   lon         the longitude matrix of your grid 
%
% OUTPUT:
%   D           the 3d matrix of images (x,y,time)
%   
%   Last modified: bgetraer@princeton.edu, 7/5/2018
%   See also: BOXGREENLAND.m, imageryseq.m, GRACE2PLMT, PLM2XYZ

% size and file name of the image sequence matrix
shp = size(lon);
filename = 'im_seqSH';

% initialize the blank 3d matrix (x, y, time)
D = zeros([shp,size(thedates,2)]);

% For each date, evaluate the spherical harmonic equations on our grid
%   this is the time consuming bit... plm2xyz takes a while
for i=1:size(thedates,2)
    % Print where we are in the loop
    fprintf('i=%d of %d \nnow starting plm2xyz\n',i,size(thedates,2))
    % Do the evaluation on the grid
    [data]=plm2xyz(squeeze(coffs(i,:,:)),lat(:),lon(:)); % vector of solutions
    D(:,:,i) = reshape(data,shp);  % put the vector into the image matrix
    % Save the images so that if MATLAB crashes you can modify your script
    %   and not repeat the completed dates.
    save(fullfile(datadir,filename),'D','thedates')
end

fprintf('PLM2GRID is complete')
end

