function [ dailyNAOvalue ] = readDailyNAO( date, smoothing )
% READDAILYNAO returns daily NAO value for requested date, with a moving
%   average if desired.
%   Opens daily NAO values from the values stored in datafiles/NAO_daily.ascii
%   and smooths with a moving average.
%
% INPUT
%   date        single or array of datenumbers
%   smoothing   integer for window of moving average (in days) (def=1)
%
% OUTPUT
%   dailyNAOvalue   single or array of NAO index values corresponding to
%                       the dates
%
% SEE ALSO:
%   MOVMEAN
%
% Last modified by bgetraer@princeton.edu 2/15/2019

% default to no smoothing
defval('smoothing',1)

dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
addpath(fullfile(dir,'/datafiles/NAO'));

% read in and smooth the NAO daily index data
NAO_daily = load('NAO_daily.ascii');
NAOdates = datenum(NAO_daily(:,1:3));
NAOdata = movmean(NAO_daily(:,4),smoothing);

% initialize output array length
dailyNAOvalue = zeros(1,length(date));

% select only the dates requested
for i = 1:length(date)
    dailyNAOvalue(i) = NAOdata(NAOdates==date(i));
end

end

