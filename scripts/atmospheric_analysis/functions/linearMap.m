function [ slopes,F ] = linearMap( data, alldates, degree )
%LINEARMAP Takes a 2d matrix timeseries, removes seasonal periodic trend,
% returns map of best fit linear model coefficient at each pixel 
%   Detailed explanation goes here
% NORMALIZE PIXEL BY PIXEL

defval('degree',1);

% The GRIS Merra nan mask
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
datadir = fullfile(dir,'datafiles/');
load(fullfile(datadir,'comparemassMerra'),'indexICEmerra','indexICE')
if size(indexICEmerra') == size(data(:,:,1))
    GRISnan = double(indexICEmerra);
    GRISnan(~indexICEmerra)=nan;
    GRISnan = GRISnan';
elseif size(indexICE) == size(data(:,:,1))
    GRISnan = double(indexICE);
    GRISnan(~indexICE)=nan;
    GRISnan = ones(size(data(:,:,1)));
else
    GRISnan = ones(size(data(:,:,1)));
    fprintf('cannot find mask of correct size, unmasked');
end

% over the ice extent
datamat = data.*GRISnan;

% vectorize it
datavec = reshape(datamat,size(datamat,1)*size(datamat,2),size(datamat,3));
meandatavec = nanmean(datavec,2);

% initialize
ml = zeros(degree+1,size(datamat,1)*size(datamat,2));

F = nan(size(datamat));
[rows,cols] = find(F);
% go through each point on the map, correct for seasonal trend, take linear
% model
warning('off','MATLAB:nearlySingularMatrix')
for i = 1:size(datamat,1)*size(datamat,2)
    if ~isnan(meandatavec(i))
        [~,fp] = periodic_m(alldates,datavec(i,:)-meandatavec(i),[365,365/2]);
        x = alldates;
        y = datavec(i,:)-meandatavec(i)-fp';
        [ml(:,i),f] = linear_m(x,y,degree);
        F(rows(i),cols(i),:) = f;
    end
end
slopes = reshape(ml(degree+1,:),size(datamat,1),size(datamat,2)).*GRISnan;
end

