function [ slopes ] = linearMap( data, alldates )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% NORMALIZE PIXEL BY PIXEL

% The GRIS Merra nan mask
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
datadir = fullfile(dir,'datafiles/');
load(fullfile(datadir,'comparemassMerra'),'indexICEmerra')
GRISnan = double(indexICEmerra);
GRISnan(~indexICEmerra)=nan;

% over the ice extent
datamat = data.*GRISnan';

% vectorize it
datavec = reshape(datamat,size(datamat,1)*size(datamat,2),size(datamat,3));
meandatavec = nanmean(datavec,2);

% initialize
ml = zeros(2,size(datamat,1)*size(datamat,2));

% go through each point on the map, correct for seasonal trend, take linear
% model
warning('off','MATLAB:nearlySingularMatrix')
for i = 1:size(datamat,1)*size(datamat,2)
    if ~isnan(meandatavec(i))
        [~,fp] = periodic_m(alldates,datavec(i,:)-meandatavec(i),[365,365/2]);
        x = alldates;
        y = datavec(i,:)-meandatavec(i)-fp';
        [ml(:,i)] = linear_m(x,y);
    end
end
slopes = reshape(ml(2,:),size(datamat,1),size(datamat,2)).*GRISnan';
end

