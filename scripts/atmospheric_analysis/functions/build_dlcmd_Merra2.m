function [ thecmd ] = build_dlcmd_Merra2(thedate)
%BUILD_DLCMD_MERRA2   Build the download command for the Merra2 data we want.
%   You give the date, we build the wget command.
%   INPUT
%       thedate     datetime 
%   OUTPUT
%       thecmd      command string
%
% Last modified by bgetraer@princeton.edu

% netCDF data file directory
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
merra_dir = fullfile(dir,'datafiles/MERRA2');
% the data url on nasa.gov
[get_url, thefilename] = build_url_Merra2(thedate);
thefiledir = fullfile(merra_dir,thefilename);
% wget download information
username = 'bgetraer';
password = 'hujguG-3mocgy-kixhiw';
cmd_format = '/usr/local/bin/wget --user %s --password %s -O %s \"%s\"';
% the system command
thecmd = sprintf(cmd_format,username,password,thefiledir,get_url);
end

