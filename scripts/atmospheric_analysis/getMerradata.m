%**************************************************************************
% Code for downloading Merra2 data netCDF files, from a text file of links.
% Merra2 Data from:
% https://disc.gsfc.nasa.gov/datasets/M2I6NPANA_V5.12.4/summary?keywords=%22MERRA-2%22
%
% execute downloads using
%**************************************************************************

days = 8:9;
filetext = '';

% write a download command for each day you want
for i=1:length(days)
    thedate = datetime(2012,7,days(i));
    cmd = build_dlcmd_Merra2(thedate);
    filetext = [filetext cmd ';\n'];
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