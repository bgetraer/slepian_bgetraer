function [ hbar ] = histBAR( hval,edges,color,sign )
%HISTBAR Creates a histogram, but formatted as a bar graph for ease of
%   plotting.
%
% INPUT
%   hval    values from histcounts
%   edges   edges from histcounts
%   color   if desired
%   sign    to flip the histogram upside down, sign
%
% OUTPUT
%   hbar    a handle to the bar graph
%
%   bgetraer@princeton.edu

defval('color','flat')
defval('sign',1)

htick = edges(1:end-1)+diff(edges)/2;
hbar = bar(htick,sign*hval,1,'FaceAlpha',0.3,'FaceColor', color);
end

