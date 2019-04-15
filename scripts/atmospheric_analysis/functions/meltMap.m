function [ filename ] = meltMap( startdate, enddate, Temp, matDir )
%MELTMAP Loads surface temperature files for the date range asked for,
% calculates the surface temperature using floorData for each day, and
% saves the matrices in a matlab data file.
%
% INPUT
%   startdate, enddate      define the date range requested, datenumbers
%   Temp                    the temperature threshold (K)
%   matDIR                  where the matlab datafiles are located
% OUTPUT
%   filename                this is the name of the new file made
% SEE ALSO:
%   FLOORDATA
%
% Last modified by bgetraer@princeton.edu 3/12/2019

datadir = fullfile('/Users/benjamingetraer/Documents/IndependentWork',...
    'slepian_bgetraer/datafiles/MERRA2');
defval('matDir', fullfile(datadir,'MerraMat'))
defval('Temp',273.15)

% the date range
thedate = datenum(startdate):1:datenum(enddate);
yrs = unique(year(thedate));

% peek at matrix size
surfTfile = sprintf('floorT%i',yrs(1));
load(fullfile(matDir,surfTfile));
% initialize the full matrix
allmeltMap = zeros(size(floorT,1),size(floorT,2),length(thedate),...
    'logical');
clear thetimedata floorT

% file name and location
if Temp == 273.15
    filename = fullfile(matDir,'allmeltMap.mat');
else
    filename = fullfile(matDir,sprintf('allmeltMap%s.mat',num2str(Temp)));
end
dateindex = 0;

if ~exist(filename,'file')
    for i = 1:length(yrs)
        y = yrs(i);
        % get the data for that year
        surfTfile = sprintf('floorT%i',y);
        load(fullfile(matDir,surfTfile))
        for j = 1:size(floorT,3)
            dateindex = dateindex + 1;
            allmeltMap(:,:,dateindex) = meltThresh(floorT(:,:,j),[],Temp);
        end
        clear thetimedata floorT
    end
    save(filename,'allmeltMap');
end

end