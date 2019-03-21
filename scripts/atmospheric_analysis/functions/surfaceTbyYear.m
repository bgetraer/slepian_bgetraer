function [ filenames, yrs, thespacelim ] = ...
    surfaceTbyYear( startdate, enddate, ncdfDir, matDir )
%SURFACETBYYEAR Loads netCDF files for the date range asked for,
% calculates the surface temperature using floorData for each day, and
% saves the matrices in a matlab data file by year.
%   This function plays "Bound for the Floor" by Local H if the song is
%   available in the matlab data file directory.
%
% INPUT
%   startdate, enddate      define the date range requested, datenumbers
%   ncdfDIR                 where the netCDF files are kept
%   matDIR                  where the new matlab datafiles should be put
% OUTPUT
%   filenames               these are the locations of the new datafiles
%                               created
%   fileyear                the year corresponding to the filename
%   thespacelim             the spatial limits of the files
% SEE ALSO:
%   GETMERRA2VAR, FLOORDATA
%
% Last modified by bgetraer@princeton.edu 3/12/2019

datadir = fullfile('/Users/benjamingetraer/Documents/IndependentWork',...
    'slepian_bgetraer/datafiles/MERRA2');
defval('ncdfDir', fullfile(datadir,'MerraGrnlandTMEAN'))
defval('matDir', fullfile(datadir,'MerraMat'))


% the date range
thedate = datenum(startdate):1:datenum(enddate);
yrs = unique(year(thedate));

% file name and location
filenames = cellstr(strcat(matDir,'/floorT',num2str(yrs(:)),'.mat'));

% get the data for the space limits once
[~,~,thespacelim] = ...
    getMerra2VAR(thedate(1),'T',ncdfDir);

if ~exist(filenames{1},'file')
    if lengthyProcessFlag('surfaceTbyYear')
        % BOUND FOR THE FLOOR
%         boundMP3 = fullfile(matDir,'Bound4theFloor.mp3');
%         if exist(boundMP3,'file')
%             [music,samplerate] = audioread(boundMP3);
%             player = audioplayer(music,samplerate);
%             play(player);
%         end
        
        for i = 1:length(yrs)
            y = yrs(i);
            if ~exist(filenames{i},'file')
                % get the data for that year
                [tdata,thetimedata] = ...
                    getMerra2VAR(thedate(year(thedate)==y),'T',ncdfDir); %#ok<ASGLU>
                % initialize the floorT matrix
                floorT = zeros(size(tdata,1),size(tdata,2),size(tdata,4));
                for j = 1:size(tdata,4)
                    floorT(:,:,j) = floorData(tdata(:,:,:,j));
                end
                
                save(filenames{i},'floorT','thetimedata')
                
                clear tdata thetimedata floorT
            end
        end
    end
end
end

