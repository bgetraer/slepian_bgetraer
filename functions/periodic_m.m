function [m,f,Rsq,power_m] = periodic_m(x,y,periods,est_error)
%PERIODIC_M is a linear model of sine and cosine terms of given periods for
%   the data y over the domain x.
%
%   Some code and logic adapted from SLEPT2RESID.m by Chris Harig


x = x(:);
y = y(:);
% some useful constants
twopi = 2*pi;
% length of data
nx=length(x);
% The frequencies being fitted in [1/time]
omega = 1./periods;
% Number of periodic components
lomega=length(omega);

% Angular frequency of the periodic terms
theta = repmat(omega,nx,1)*twopi;   % the angular frequency terms
exes = repmat(x,1,lomega);          % the domain terms
th_x = theta.*exes;                 % the full terms

% Here is our linear operator: first all the cosine terms then the sine
%   terms
G = [cos(th_x) sin(th_x)];

% Are there a-priori error estimates?
if exist('est_error','var')
    % Create a weighting matrix, weight the G and x matrices.
    W = diag(1./est_error);
    G = W*G1;
    y = W*y;
end

% Do the fitting by minimizing least squares. Unweighted, or weighted by
% the estimated errors
m = (G'*G)\G'*y;   % the model parameters
f = G*m;            % the modeled signal

% return variance reduction
Rsq = 1-var(y-f)/var(y);

% return the combined power of sin & cos for each period
stacked_m = [m(1:end/2)';m(end/2+1:end)']; %stack corresponding coeffs.
power_m = sum(stacked_m.^2,1);  %square and sum corresponding coeffs.
end