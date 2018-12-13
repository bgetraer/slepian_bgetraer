function [ range ] = ncdf_timerange(filename)
date1 = ncreadatt(filename,'/','RangeBeginningDate'); 
time1 = ncreadatt(filename,'/','RangeBeginningTime'); 
date2 = ncreadatt(filename,'/','RangeEndingDate'); 
time2 = ncreadatt(filename,'/','RangeEndingTime'); 

range{1} = datestr(datenum(horzcat(date1,',',time1)),'dd-mmm-yyyy,HH:MM');
range{2} = datestr(datenum(horzcat(date2,',',time2)),'dd-mmm-yyyy,HH:MM');
end

