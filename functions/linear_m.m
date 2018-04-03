function [m,f,DELTA,mint,t] = linear_m(x,y,n)
%LINEAR_M returns an unweighted linear model of a 1st or 2nd degree 
%   polynomial of form y = m1 + m2*x + 1/2*m3*x^2. Operates like POLYFIT, 
%   but adjusts the parameter output for m3 scaled by 1/2 to represent
%   acceleration, and returns the modeled signal f(t). if x is longer than
%   y, assumes that extra x terms are appended on the end for extrapolation
%   of the fit.
%
%
%INPUTS
%   x   the independent variable
%   y   the dependent variable
%   n   degree (1 or 2) default: 1
%
%OUTPUTS
%   m   the model parameters
%   f   the modeled dependent variable
%   DELTA   the estimated predictive range of the regression 
%   t   the independent variable (same as x)
%
% See also: FIT POLYFIT POLYVAL
%
% Last modified by bgetraer@princeton.edu 1/8/2018

% put them in column form
xprime = x(:);
y = y(:);
x = xprime(1:length(y));

% defval('n',1)
if n>2
    error('Only supports degree 1 and 2 polynomial regression');
end

% Model: 1st or 2nd degree polynomial regression
%   y(t) = m(1) + m(2)*x + 1/2*m(3)*x^2 
%   y = G*m for data y, linear operator G, and model parameters m

G = ones(length(x),n+1);      % constant terms
G(:,2) = x;                 % linear terms
if n==2
        G(:,3) = 1/2.*x.^2;         % quadratic terms
end


%OUTPUTS
% model parameters
m = (G'*G)\G'*y;

%GET ESTIMATES ON ERROR
% note: P should be the same as m but with a factor of 1/2 in the m(3)
% and rescaled for better error estimation
[P,S,mu] = polyfit(x,y,n);
[f,DELTA] = polyval(P,xprime,S,mu);

%GET ESTIMATES ON CONFIDENCE INTERVAL OF COEFFICIENTS
% note: fitresult should return a model equivalent to m but with a factor 
% of 1/2 in the m(3)

fitresult = fit(x,y,sprintf('poly%i',n));
ci = confint(fitresult);
mint = zeros(size(m));
for i = 1:n+1
    mint(i) = abs(ci(1,end+1-i)-m(i));
end

% independent variable
t = xprime;
end

