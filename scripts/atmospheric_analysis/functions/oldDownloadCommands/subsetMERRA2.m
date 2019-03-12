function [ subset ] = subsetMERRA2(varname, time, lat, lon, lev )
%SUBSETMERRA2 Builds the subset arguments for the Merra2 download link

if varname == 'all'
    varname = {'PS';'V';'T';'SLP';'U';'QV';'H';'O3'};
end

defval('time',[0 3]);
defval('lat',[240 360]);
defval('lon',[144 288]);
defval('lev',[0 41]);

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


% subsetinfo4 = ['.nc?'...
%     'PS[0:3][240:360][144:288],'...
%     'V[0:3][0:41][240:360][144:288],'...
%     'T[0:3][0:41][240:360][144:288],'...
%     'SLP[0:3][240:360][144:288],'...
%     'U[0:3][0:41][240:360][144:288],'...
%     'QV[0:3][0:41][240:360][144:288],'...
%     'H[0:3][0:41][240:360][144:288],'...
%     'O3[0:3][0:41][240:360][144:288],'...
%     'time[0:3],lat[240:360],lon[144:288],'...
%     'lev[0:41]'];

end

