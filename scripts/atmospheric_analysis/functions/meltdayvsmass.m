function [ thismelt, thismass ] = meltdayvsmass(mnth,yr,massimages,massdates,massX,massY,...
    meltimages,meltdates,meltX,meltY)
%MELTDAYVSMASS takes in date requests and data for MERRA and GRACE,
% returns the change in data across those dates and interpolates the MERRA
% data onto the GRACE points.
%   
% Last Modified: bgetraer@princeton.edu 3/26/2019

if length(mnth)~=2
    error('choose two months for the comparison')
elseif yr(1)==yr(end) && mnth(1)>=mnth(end)
    error('choose two different months in ascending order')
end

d1 = monthnum(mnth(1),yr(1),massdates);
d2 = monthnum(mnth(end),yr(end),massdates);
thismass = massimages(:,:,d2) - ...
    massimages(:,:,d1);

d1 = find(meltdates==datenum(yr(1),mnth(1),1));
d2 = find(meltdates==datenum(yr(end),mnth(end),1));
thismelt = imgaussfilt(sum(meltimages(:,:,d1:d2),3),1);

thismelt = interp2(meltX,meltY,thismelt',massX,massY);
end

