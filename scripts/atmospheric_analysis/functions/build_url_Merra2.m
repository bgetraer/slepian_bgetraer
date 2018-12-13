function [ theurl, thefilename] = build_url_Merra2(thedate)
%BUILD_URL_MERRA2   Build the download url of the Merra2 data we want.
%   You give the date, we build the url. In our case, we want the date
%   associated with the inst6_3d_ana_Np files.
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
url1 = 'https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I6NPANA';
datestr1 = datestr(thedate,'yyyy/mm/');
url2 = ['.5.12.4/' datestr1 'MERRA2_400.'];
datestr2 = datestr(thedate,'yyyymmdd');
filename3 = sprintf('inst6_3d_ana_Np.%s.nc4',datestr2);
% specifies the subset of data I want
subsetinfo4 = '.nc?PS[0:3][240:360][144:288],V[0:3][0:41][240:360][144:288],T[0:3][0:41][240:360][144:288],SLP[0:3][240:360][144:288],U[0:3][0:41][240:360][144:288],QV[0:3][0:41][240:360][144:288],H[0:3][0:41][240:360][144:288],O3[0:3][0:41][240:360][144:288],time,lat[240:360],lon[144:288],lev';
theurl = [url1 url2 filename3 subsetinfo4];
thefilename = sprintf('inst6_3d_ana_Np_SUBSET.%s.nc4',datestr2);
end

