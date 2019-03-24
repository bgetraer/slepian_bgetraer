function [ filenames, yrs, thespacelim ] = ...
    surfaceTbyYear( startdate, enddate, ncdfDir, matDir )
%SURFACETBYYEAR Loads netCDF files for the date range asked for,
% calculates the surface temperature using floorData and pulls the maximum
% temperature for each day, and saves the matrices in a matlab data file by
% year.
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
                
                maxfloorT = zeros(size(tdata,1),size(tdata,2),length(thedates));
                % HERE WE NEED TO GO THROUGH EACH DAY, AND FOR EACH DAY
                % FIND THE MAX T. Assume uniform time points per day.
                tpd = size(floorT,3)/length(thedates);
                if tpd ~= 1
                    if tpd == round(tpd)
                        for k = 1:length(thedates)
                            thisrange = (k-1)*tpd+1 : tpd+(k-1)*tpd;
                            maxfloorT(:,:,k) = max(floorT(:,:,thisrange),[],3);
                        end
                    else
                        error('non-uniform time points per day')
                    end
                else
                    maxfloorT = floorT;
                end
                save(filenames{i},'maxfloorT','thetimedata')
                
                clear tdata thetimedata floorT
            end
        end
    end
end
end

