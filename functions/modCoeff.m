function [sf,vrnc] = modCoeff(thedates,y)
%MODCOEFF Models coefficients, here designed for my wavelet coefficients,
%using constant, linear, quadratic, annual and semiannual terms. We do an
%unweighted inversion, and return the modeled signal and the variance of the
%errors. Inspired by the modeling of Slepian coefficients.
%
% bgetraer@princeton.edu 5/1/2018

recordlength = (thedates(end)-thedates(1));

x = (thedates-thedates(1))/recordlength; % scale 0 to 1

yearlength = 365.2422; % in days
freqANNUAL = yearlength/recordlength;
freqSEMI = yearlength/2/recordlength;

freq8 = yearlength*8/recordlength;


omega = 2*pi/freqANNUAL;
phi = 2*pi/freqSEMI;
theta = 2*pi/freq8;

Gc = ones(1,length(x));
Glin = x;
Gquad = x.^2;
Gcub = x.^3;
GsinA = sin(omega.*x);
GcosA = cos(omega.*x);
GsinS = sin(phi.*x);
GcosS = cos(phi.*x);
Gsin8 = sin(theta.*x);
Gcos8 = cos(theta.*x);

G = [Gc(:) Glin(:) Gquad(:) Gcub(:) GsinA(:) GcosA(:) GsinS(:) GcosS(:)];

% G = [Gc(:) Glin(:) Gquad(:) Gcub(:) GsinA(:) GcosA(:) GsinS(:) GcosS(:) Gsin8(:) Gcos8(:)];


mod = G\y;

sf = G*mod;

vrnc = var(y-sf);



end

