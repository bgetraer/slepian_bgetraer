function shpower(l,m,plot)
%SHPOWER plots the spatial distribution of power for a given Spherical
%Harmonic of order l, degree m
%
%INPUT
%   l       a single angular degree (default=random between 1:60)
%   m       a single angular order (default=random between 0:l)
%   plot    
%           1   Mollweide projection
%           2   3-D sphere, no topography
%           3   Flat rectangular projection for entire globe
%
%Last modified by bgetraer@princeton.edu, 12/03/2017

defval('l',randi(60));
defval('m',randi(l));

if l<1
    error('l must be a single angular degree greater than 0')
end

if m>l
    error('m must be a single angular order less than or equal to l')
end

[~,~,~,blank]=addmon(l);
blank(blank(:,1)==l&blank(:,2)==m,3:4) = [1,1];


switch plot
    case 1
        plotplm(blank,[],[],7,1)
    case 2
        plotplm(blank,[],[],8,1)
    case 3
        plotplm(blank,[],[],9,1)
    otherwise
        error('Not a valid plotting method')
end
title(strcat('$l=$',sprintf('%i',l),' $m=$',sprintf('%i',m)),'interpreter','latex');
end

