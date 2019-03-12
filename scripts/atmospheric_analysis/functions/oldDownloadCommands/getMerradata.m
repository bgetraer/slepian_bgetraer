%**************************************************************************
% Process for downloading Merra2 data netCDF files, from a text file of links.
% Merra2 Data from:
% https://disc.gsfc.nasa.gov/datasets/M2I6NPANA_V5.12.4/summary?keywords=%22MERRA-2%22
%
% execute downloads using
%**************************************************************************

% example netCDF file directory
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
merradata_dir = fullfile(dir,'datafiles/MERRA2');
addpath('functions')
addpath(fullfile(dir,'functions'))
setworkspace('/Users/benjamingetraer/Documents/IndependentWork/SH_Workspace');

%%
date1 = datetime(2003,2,1);
date2 = datetime(2010,2,1);
daynum = date2;%:date2;
filetext = '# MERRA2 DOWNLOAD REQUEST \n';


% write a download command for each day you want
for i=1:length(daynum)
    thedate = daynum(i);
    thefilename = sprintf(...
        'inst6_3d_ana_Np_SUBSET.%s.nc4',datestr(thedate,'yyyymmdd'));
    if ~exist(fullfile(merradata_dir,thefilename),'file')
        
        cmd = build_dlcmd_Merra2(thedate,{'T'});
        filetext = [filetext cmd ';\n'];
    end
end
% create text file
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
merra_dir = fullfile(dir,'scripts/atmospheric_analysis');
filelocation = fullfile(merra_dir,'dlMerra2.txt');
fid = fopen(filelocation,'wt');
fprintf(fid, filetext);
fclose(fid);
%copy shell script command
clipboard('copy', ['sh -e ' filelocation])

%%
wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies --content-disposition
