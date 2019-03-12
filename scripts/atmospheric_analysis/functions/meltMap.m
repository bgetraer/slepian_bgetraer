function [ filename ] = meltMap( startdate, enddate, matDir )
%MELTMAP Loads surface temperature files for the date range asked for,
% calculates the surface temperature using floorData for each day, and
% saves the matrices in a matlab data file.
%
% INPUT
%   startdate, enddate      define the date range requested, datenumbers
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
filename = fullfile(matDir,'allmeltMap.mat');
dateindex = 0;

if ~exist(filename,'file')
    for i = 1:length(yrs)
        y = yrs(i);
        % get the data for that year
        surfTfile = sprintf('floorT%i',y);
        load(fullfile(matDir,surfTfile))
        for j = 1:size(floorT,3)
            dateindex = dateindex + 1;
            allmeltMap(:,:,dateindex) = meltThresh(floorT(:,:,j));
        end
        clear thetimedata floorT
    end
    save(filename,'allmeltMap');
end

end


%% OLD CODE, USED TO MAKE FILES BY YEAR
% filenames = cell(1,length(yrs));
% fileyear = yrs;
% for i = 1:length(yrs)
%     y = yrs(i);
%     % file name and location
%     fileName = sprintf('meltMap%i',y);
%     filenames{i} = fullfile(matDir,fileName);
%     
%     if ~exist(filenames{i},'file')
%         % get the data for that year
%         surfTfile = sprintf('floorT%i',y);
%         load(fullfile(matDir,surfTfile))
%         % initialize the meltmap matrix
%         meltmap = zeros(size(floorT,1),size(floorT,2),size(floorT,3));
%         for j = 1:size(floorT,3)
%             meltmap(:,:,j) = meltThresh(floorT(:,:,j));
%             
%         end
%         save(fullfile(matDir,fileName),'meltmap','thetimedata')
%         
%         clear meltmap thetimedata floorT
%     end
% end


