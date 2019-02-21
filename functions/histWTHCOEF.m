function [recon_image, prctGone, nLeft] = histWTHCOEF( originalimage, target, wname )
%HISTWTHCOEF Plots a histogram of wavelet coefficients, with a visual
% representation of the thresholding done in PRCTILETHOLD.
%   Displays a histogram of the wavelet coefficients by magnitude,
%   distinguishing between those below and above the percentile threshold.
%
% INPUT
%   originalimage   the image which is being deconstructed
%   target          the value of bias or invar to be maintained (def. 90%)
%   wname           wavelet to be used in the decomposition (def. haar)
%
% OUTPUT
%   recon_image     the reconstructed coefficient-thresholded image
%   prctGone        the percentage of coefficients eliminated
%   nLeft           the number of coefficients remaining
%
% SEE ALSO:
%   PRCTILETHOLD
%
% Last modified: bgetraer@princeton.edu 2/21/2019

defval('target',0.9)
defval('wname','haar')

% Calculate threshold and reconstruct image
[ ~, ~, ~, CT, S , C, T] = prctileThold(originalimage, target,wname); 
recon_image = waverec2(CT,S,wname);

% Plot
figure(1); cla
% order the coefficients least to greatest and take the log
Cdesc = sort(abs(C),'descend');
X = log(Cdesc);

hold on
% Histogram plot:
[y, n] = hist(X,100); % y: values; n: bin centers
yn = y/sum(y);
lT = log(T);
ind = n>lT; % bin centers: greater or smaller than threshold?
h_more = bar(n(ind), yn(ind), 1, 'b'); % for greater: use blue
h_less = bar(n(~ind), yn(~ind), 1, 'r','edgecolor','none'); % for smaller: use red

% locate bar around threshold: it needs the two colors
[~, nd] = min(abs(n-lT)); 
patch([(n(nd+1)+n(nd))/2 lT lT (n(nd+1)+n(nd))/2], [0 0 yn(nd) yn(nd)], 'b');

% draw line at threshold
h_line = plot(log([T T]),ylim,'k');

% explain the plot a bit
xlabel('$log(\big|$wavelet coefficients$\big|)$','interpreter','latex')
ylabel('Frequency','interpreter','latex')
lgd = legend([h_less, h_more, h_line],'Below threshold', 'Above threshold',...
    sprintf('Threshold $=e^{%0.3f}$',lT));
set(lgd,'interpreter','latex','location','northwest')
set(gca,'fontsize',12)

% percentage of coefficients thrown out:
prctGone = sum(yn(~ind));
% number left
nLeft = sum(CT(:)~=0);
end

