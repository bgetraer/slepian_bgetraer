%%*************************************************************************
% CIRCULAR SLEPIAN BASES 
% Circular bases are designated about specific points of interest around the
% globe - major oceans, continents, etc. - and slepian functions are
% calculated for those areas. 
%% SETUP WORKSPACE
addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

radius = 15;

% plot location of circular bases over major oceans and continents
figure(1)
subplot(2,1,1)
setcircbases(1,radius,1);
subplot(2,1,2)
setcircbases(2,radius,1);

%% TEST TRANSFORMATION TO SLEPIAN BASIS
centertype = 1;
% define the centers, the bandlimit, and the radius
centers = setcircbases(centertype,radius);
bandlimit = 20;
% define the SH coefficients corresponding to the Slepian basis
[G, N] = slepcircbases(radius,bandlimit,centers);

% Take some of the Slepian functions from each basis plot them
figure(2);clf;hold on;
[lmcosi_mat,lmcosi_sum,alphas] = grabanalpha(G,N,bandlimit,1,3);

setcircbases(centertype,radius,1,0);
for i=1:length(centers)
    text(centers(1,i),centers(2,i),strcat('$\alpha=$',sprintf('%i',alphas(i))),...
        'HorizontalAlignment','center','interpreter','latex');
end

