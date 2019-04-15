function [month_num, exact_date_num] = monthnum(MONTH,YEAR,thedates)
%MONTHNUM This function takes the date vector from grace2plmt and  
% returns the index corresponding to a requested month and year, or the 
% closest month on record.           
%   MONTH is a single month requested, as an integer 1:12
%   YEAR is a single year requested, as an integer 2002:2017
%
%   month_num is a single integer that is the index corresponding to a 
%       requested month and year.
%   exact_date_num is the MATLAB datenumber corresponding to the median day
%       of the given monthly measurement.
%
% example:
%       [potcoffs,cal_errors,thedates] = grace2plmt('CSR','RL05','SD',0);
%       [month_num, exact_date_num] = plmt_monthnum(MONTH,YEAR,thedates);
%       MONTH_YEAR_COEFF = potcoffs(month_num,:,:)
%
% SEE ALSO:
%   DATENUM, 
%
% 12.28.2017 bgetraer@princeton.edu

% ensure valid date request
if ~any(MONTH==1:12),error('MONTH is not an integer 1:12');end
if ~any(YEAR==2002:2017),...
        error('YEAR is not an integer within years of record');end
if ~any(month(thedates)==MONTH & year(thedates)==YEAR),...
        warning('MONTH, YEAR is not within provided record');end
if any(month(thedates)==MONTH & year(thedates)==YEAR)
    the_date_num = thedates(month(thedates)==MONTH & year(thedates)==YEAR);
else, the_date_num = datenum(YEAR,MONTH,15);
end

% find and return nearest date
[~, month_num] = min(abs(thedates-the_date_num(1)));
exact_date_num = thedates(month_num);
end

