addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
addpath(datadir)
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

load('Greenland60data');
% find peaks in Greenland signal
figure(1)
clf
plot(thedates,ESTtotal)
findpeaks(ESTtotal)
% store index
[pk,loc]=findpeaks(ESTtotal);

peakdates = thedates(loc);

month(peakdates)

% how mass has changed
dp = diff(pk);
mean(dp)
std(dp)
