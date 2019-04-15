function [ isMeltingMap ] = meltThresh( tdata,filter,meltThreshold )
%MELTTHRESH Returns a binary map thresholded at the freezing point.
%
%INPUT
%   tdata   temperature data matrix. If 3d, floorData will take the first 
%               real values of the stack of layers and filter them to
%               remove artifacts of the wierd stepping between pressure
%               levels.
%OUTPUT
%   isMeltingMap    a binary map thresholded at the freezing point.
%
% Last modified by bgetraer@princeton.edu 3/12/2019

defval('filter',1)
defval('meltThreshold',273.16)


if length(size(tdata)) == 3
    t_floor = floorData(tdata);
    t_filter = imgaussfilt(t_floor,filter);
    isMeltingMap = t_filter > meltThreshold;
else
    isMeltingMap = tdata > meltThreshold;
end
end

