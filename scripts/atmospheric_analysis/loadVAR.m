%**************************************************************************
% Script for loading data over the Greenland Ice Sheet
% Merra2 Data from:
% https://disc.gsfc.nasa.gov/datasets/
%
% Last modified by: bgetraer@princeton.edu 3/30/2019
%**************************************************************************

%% locate slepian_bgetraer function and datafile directories, and set workspace
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
merra_dir = fullfile(dir,'datafiles/MERRA2');
addpath(fullfile(dir,'scripts/atmospheric_analysis/functions'))
addpath(fullfile(dir,'functions'))
setworkspace('/Users/benjamingetraer/Documents/IndependentWork/SH_Workspace');
matDir = fullfile(merra_dir,'MerraMat');

%%
startdate = datenum(2003,01,01);
enddate = datenum(2017,12,31);
alldates = startdate:enddate;

% dataset = 'tavg1_2d_slv_Nx';
% dataset = 'instM_2d_asm_Nx';
% dataset = 'tavg1_2d_rad_Nx';
% dataset = 'statM_2d_slv_Nx';
% dir = fullfile(merra_dir, 'MerraRadiation');
dir = fullfile(merra_dir, 'MerraGrnlandTMEAN');
variable = '';
matfilename = fullfile(matDir,sprintf('%s2003-2017',variable));

if lengthyProcessFlag
[ data, t, spacelim ] = ...
    getMerra2VAR( alldates, variable, dir ,dataset);
squeeze(data);
save(matfilename,'data','t','spacelim')
end

