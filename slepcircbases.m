function [ G,N ] = slepcircbases( r,L,centers )
%SLEPCIRCBASES 
% Construct slepian expansion GLMALPHAPTO for rotated cap centers of given
% radius and bandlimit. Outputs an (lm)X(alpha)X(basis) matrix with
% spherical harmonic coefficients of the Slepian functions, where lm are 
% coefficients, alpha are the eigentapers, and basis are the circular caps 
% rotated to a given center.
% 
%INPUT
%   r           Angular extent of the spherical cap (degrees)
%   L           Bandwidth (maximum angular degree)
%   centers     array of geographic centers in [lon;lat]
%                   or
%               1   ocean centers
%               2   continent centers
%                       default = 1
%** NOTE **
%centers must have latitude (-90 to 90) NOT co-latitude (0 to 180)
%
%OUTPUT
%   G           The (lm)X(alpha)X(basis) matrix of coefficients
%   N           The vector of corresponding shannon numbers
%       See also:
%           GLMALPHAPTO SETCIRCBASES CAPLOC
%
%Last modified by bgetraer@princeton.edu, 11/24/2017

defval('centers',1);
if centers == 1
        centers = setcircbases(1);
elseif centers == 2
        centers = setcircbases(2);
end

% Create blank G
[~,~,~,~,~,~,~,~,~,ronm]=addmon(L);
G = zeros([length(ronm),length(ronm),length(centers)]);

% Find the Slepian Basis for centers
for j = 1:size(centers,2)
    centers(2,j)
    % get Linear Combination of SH functions for each rank Slepian Function
    [Gcap,V,EL,EM,N,GM2AL,MTAP,IMTAP] = glmalphapto(r,L,centers(1,j),90-centers(2,j));
    % save the SH combinations for each rank Slepian Function
    G(:,:,j) = Gcap;
end
end