%%*************************************************************************
% CIRCULAR SLEPIAN BASES 
% Circular bases are designated about specific points of interest around the
% globe - major oceans, continents, etc. - and slepian functions are
% calculated for those areas. 
%% SETUP WORKSPACE
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
radius = 15;
% define the SH coefficients corresponding to the Slepian basis
[G, N] = slepcircbases(radius,bandlimit,centers);
% Take some of the Slepian functions from each basis and put them into a 
% format that PLM2XYZ knows how to interpret
[lmcosi_mat,lmcosi_sum,alphas] = grabanalpha(G,N,bandlimit);

% plot
figure(2);clf;hold on;
plotplm(lmcosi_sum,[],[],7,1);
setcircbases(centertype,radius,1,0);
for i=1:length(centers)
    text(centers(1,i),centers(2,i),strcat('$\alpha=$',sprintf('%i',alphas(i))),...
        'HorizontalAlignment','center','interpreter','latex');
end

%% Set up spherical harmonic coefficients
% % import the GRACE spherical harmonic coefficients
% [potcoffs,cal_errors,thedates]=grace2plmt('CSR','RL05','SD',0);
% % isolate a single month to test
% theelsems=squeeze(potcoffs(1,:,1:2));   % preserve the lm section
% Oct = squeeze(potcoffs(plmt_monthnum(10,2004,thedates),:,3:4));
% Apr = squeeze(potcoffs(plmt_monthnum(4,2005,thedates),:,3:4));
% plm_test = [theelsems Oct];
% 
% location_index = 1;
% % transform GRACE harmonics into Slepian expansion coefficients
% [falpha,V,N,MTAP,C] = plm2slep(plm_test,radius,bandlimit,...
%         cont_centers(1,location_index),cont_centers(2,location_index));
% % Slepian expansion coefficients to the new basis
% Gfalpha = cont_G(:,:,location_index)'.*falpha;
% 
% 
% %%
% for i = 1:length(cont_centers)
%     [falpha,V,N,MTAP,C] = plm2slep(plm_test,radius,bandlimit,...
%         cont_centers(1,i),cont_centers(2,i));
%     figure(i)
%     plotplm(falpha.*C)
% end
% %%
% figure(3)
% subplot(1,2,1)
% plotplm(C)
% colorbar
% caxis([-2 2]*1e-1)
% 
% subplot(1,2,2)
% plotplm(G)
% colorbar
% caxis([-2 2]*1e-1)
% %% plot of circular bases in lat/lon cartesian space
% figure(1)
% plotcont;   % plot the continent outlines
% hold on;
% 
% [circLON circLAT] = caploc([316,72],radius,100,1);
% plot(circLON, circLAT,'.k','markersize',0.25);
