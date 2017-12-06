function [ output_args ] = sleppower(G,alpha)
%SLEPPOWER plots the spatial distribution of power for a Slepian expansion
%given by G, of order alpha
%
%INPUT
%   G       The (lm)X(alpha)X(basis) square matrix of coefficients
%   alpha   a single eigentaper order (default=random between 0:10)
%   plot    
%           1   Mollweide projection
%           2   Flat rectangular projection for entire globe
%           3   3-D sphere, no topography
%

%Last modified by bgetraer@princeton.edu, 12/3/2017

defval('alpha',randi(10));

[~,~,~,blank]=addmon(l);
blank(blank(:,1)==l&blank(:,2)==m,3:4) = [1,1];

plotplm(blank,[],[],8,1)
title(strcat('$l=$',sprintf('%i',l),' $m=$',sprintf('%i',m)),'interpreter','latex');

end

