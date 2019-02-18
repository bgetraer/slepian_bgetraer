function [ x,y ] = selectNAO( NAOdates,NAOdata,index_season,index_outlier, index_grace)
%SELECTNAO returns a selection of NAO data based on requested indices for
%season, outliers, and time-period.
% INPUT
%       index_season    =   'summer'    (may-oct)
%                           'winter'    (nov-apr)
%                           ''          (all months)
%       index_outlier   =   'outlier%fX' 
%                               %f  = factor of std(NAOdata(theIndex))
%                               X   = 'A' all, 'P' positive, 'N' negative
%                       for example: 'outlier2.5N' returns negative
%                       outliers 2.5x the std(NAOdata(theIndex)) from the
%                       average.
%       index_grace     =   'grace'     (2002 onwards)
%                           ''          (entire dataset)
%
% Last modified by bgetraer@princeton.edu


% Indexes to call for picking apart the NAO data
theIndex = NAOdata==NAOdata; % initialize index

% temporal indices
switch index_season
    case 'summer'
        theIndex = theIndex & (month(NAOdates)>=5 & month(NAOdates)<=9);
    case 'winter'
        theIndex = theIndex & (month(NAOdates)>10 | month(NAOdates)<5);
end

% value indices
if length(index_outlier)>8
    c = str2double(index_outlier(8:end-1));
    sign = index_outlier(end);
    mu = mean(NAOdata(theIndex));
    sigma = std(NAOdata(theIndex));
    switch sign
        case 'A'
            a = 1; b = 1;
        case 'P'
            a = 1; b = 0;
        case 'N'
            a = 0; b = 1;
    end
    theIndex = theIndex & ...
        (a*(NAOdata >= mu + c*sigma) | b*(NAOdata <= mu - c*sigma));
end

switch index_grace
    case 'grace'
    theIndex = theIndex & year(NAOdates)>=2002;
end


x = NAOdates(theIndex);
y = NAOdata(theIndex);

end

