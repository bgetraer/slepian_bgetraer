function [ output_args ] = periodic_m(x,y,periods )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%      periods is the period in units of x of a function

% How many data?
defval('degree',2)
defval('givenerrors',ones(size(x)));

nx=length(x);

% RESCALING to improve the solution
mu1 = mean(x); % mean
mu2 = std(x); % standard deviation
% Make a new x-vector with this information
xprime = (x - mu1)/mu2;

% The frequencies being fitted in [1/days]
omega = 1./periods(:);
% Rescale these to our new xprime
omega = omega*mu2;
% How many periodic components?
lomega=length(omega);

% We will have the same number of G matrices as order of polynomial fit.
% These matrices are smallish, so make all 3 regardless of whether you want
% them all.
G1 = []; % For line fits
G2 = []; % For quadratic fits
G3 = []; % For cubic fits
% Mean term
if degree(1) >= 0
  G1 = [G1 ones(size(xprime'))];
  G2 = [G2 ones(size(xprime'))];
  G3 = [G3 ones(size(xprime'))];
end
% Secular term
if degree(1) >= 1
  G1 = [G1 (xprime)'];
  G2 = [G2 (xprime)'];
  G3 = [G3 (xprime)'];
end
% Quadratic term
if degree(1) >= 2
  G2 = [G2 (xprime)'.^2];
  G3 = [G3 (xprime)'.^2];
end
% Cubic term
if fitwhat(1) == 3
  G3 = [G3 (xprime)'.^3];
end
% Angular frequency in radians/(rescaled day) of the periodic terms
  th_o= repmat(omega,nx,1)*2*pi.*repmat((xprime)',1,lomega);
  G1 = [G1 cos(th_o) sin(th_o)];
  G2 = [G2 cos(th_o) sin(th_o)];
  G3 = [G3 cos(th_o) sin(th_o)];
  
  %Solving
  % If we have a priori error information, create a weighting matrix, and
    % change the G and d matrices to reflect this.  Since each coefficient
    % has its own weighting, we have to invert them separately.
    W = diag([1./givenerrors(:,index)]);
    d = slept(:,index);
    G1w = W*G1;
    G2w = W*G2;
    G3w = W*G3;
    dw = W*d;
    
    lomega = length(omega);
    myomega=omega;
    th=th_o;
%%%
    % First order polynomial
    %%%
    
    % Do the fitting by minimizing least squares
    % First, the linear fit with periodics
    mL2_1 = (G1w'*G1w)\(G1w'*dw) ;
end

