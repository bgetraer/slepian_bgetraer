function [ thecmd ] = build_dlcmd_Merra2(thedate,varname, time, lat, lon, lev)
%BUILD_DLCMD_MERRA2   Build the download command for the Merra2 data we want.
%   You give the date, we build the wget command.
%   INPUT
%       thedate     datetime
%       subset
%   OUTPUT
%       thecmd      command string
%
% Last modified by bgetraer@princeton.edu 3/1/2019

if varname{1} == 'all'
    varname = {'PS';'V';'T';'SLP';'U';'QV';'H';'O3'};
end

defval('time',[0 3]);
defval('lat',[240 360]);
defval('lon',[144 288]);
defval('lev',[0 41]);

% netCDF data file directory
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
merra_dir = fullfile(dir,'datafiles/MERRA2');
% the data url on nasa.gov
[get_url, thefilename] = build_url_Merra2(thedate,varname, time, lat, lon, lev);
thefiledir = fullfile(merra_dir,thefilename);
% wget download information
username = 'bgetraer';
password = 'hujguG-3mocgy-kixhiw';
cmd_format = '/usr/local/bin/wget --user %s --password %s -O %s \"%s\"';
% the system command
thecmd = sprintf(cmd_format,username,password,thefiledir,get_url);
end

% *************************************************************************
% *************************************************************************
function [ theurl, thefilename] = ...
    build_url_Merra2(thedate,varname, time, lat, lon, lev)
%BUILD_URL_MERRA2   Build the download url of the Merra2 data we want.
%   You give the date, we build the url. In our case, we want the date
%   associated with the inst6_3d_ana_Np files.
%
% See:
% https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I6NPANA.5.12.4/
%
%   URL format:
%       'https://DIRECTORY/FILETYPE.YYYYMMDD.nc4.nc?SUBSETINFO[]'
%   INPUT
%       thedate     datetime
%   OUTPUT
%       theurl      url string
%
% Last modified by bgetraer@princeton.edu

% the data url on nasa.gov
url1 = ['https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/MERRA2/',...
    'M2I6NPANA.5.12.4/'];
datestr1 = datestr(thedate,'yyyy/mm/');
yr = year(thedate);
if yr < 1992
        version = 'MERRA2_100.';
else if yr < 2001
        version = 'MERRA2_200.';
    else if yr < 2011
        version = 'MERRA2_300.';
        else if yr >= 2011
        version = 'MERRA2_400.';
            end
        end
    end
end
url2 = [datestr1 version];
datestr2 = datestr(thedate,'yyyymmdd');
filename3 = sprintf('inst6_3d_ana_Np.%s.nc4',datestr2);
% specifies the subset of data I want

subsetinfo4 = subsetMERRA2(varname, time, lat, lon, lev);

theurl = [url1 url2 filename3 subsetinfo4];
thefilename = sprintf('inst6_3d_ana_Np_SUBSET.%s.nc4',datestr2);
end

% *************************************************************************
% *************************************************************************
function [ subset ] = subsetMERRA2(varname, time, lat, lon, lev )
%SUBSETMERRA2 Builds the subset arguments for the Merra2 download link

formatString = '[%i:%i]';
tString = sprintf(formatString,time);
latString = sprintf(formatString,lat);
lonString = sprintf(formatString,lon);
levString = sprintf(formatString,lev);

subsetString1 = [tString lonString latString];
subsetString2 = [tString levString lonString latString];

varString = cell(1,length(varname));

subset = '.nc?';

for i = 1:length(varname)
    if strcmp(varname{i},'PS') || strcmp(varname{i}, 'SLP')
        varString{i} = strcat(varname{i}, subsetString1);
    else
        varString{i} = strcat(varname{i}, subsetString2);
    end
    subset = strcat(subset, varString{i},',');
end

subset = strcat(subset, 'time', tString, ',lat',latString,...
    ',lon', lonString, ',lev',levString);
end



