function [sigma,mu,stdline_handle,meanline_handle] = stdline(x,y,c,style)
% STDLINE horizontal dashed line illustrating one std of y in the color c
%   and linestyle style around the mean
%
% INPUTS
%   x,y         where is the data
%   c,style     how do you want the line formatted
%
% OUTPUTS
%   sigma,mu        std() and mean()
%   stdline_handle  line handle for legend
%
% SEE ALSO: 
%   STD, MEAN, LINE
%
% Last modified by bgetraer@princeton.edu

defval('c','k')
defval('style','--')
mu = mean(y);
sigma = std(y);
hold on;
line([x(1),x(end)],[mu-sigma,mu-sigma],'color',c,'linestyle',style);
stdline_handle = line([x(1),x(end)],[mu+sigma,mu+sigma],'color',c,'linestyle',style);
% meanline_handle = line([x(1),x(end)],[mu,mu],'color',c,'linestyle','--');
end