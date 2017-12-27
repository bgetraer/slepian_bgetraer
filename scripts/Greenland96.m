%**************************************************************************
% GREENLAND 96 - expand spherical harmonics of bandwidth L=96 into the
% bandwidth L=90 slepian basis for Greenland.
%
% NOTE: this script makes use of the functions GRACE2SLEPT_INPROGRESS and
% GRACE2PLMT_INPROGRESS which are modified to handle the CSR RL05 96X96
% data product.
%
% Last modified by bgetraer@princeton.edu, 12/25/2017
%**************************************************************************
addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');

% load the spherical harmonic coefficients
Ldata = 96; % bandwidth of the GRACE data
[dcoffs,sh_calerrors,~]=grace2plmt_inprogress('CSR','RL05','SD',0,Ldata);
%% expand into the first N eigentapers of the Greenland 90x90 basis
L = 90; % bandwidth of the basis into which we want to expand
J = 'N';    % how many eigentapers we want
[slepcoffs,slep_calerrors,thedates,TH,G,CC,V,N] = ...
    grace2slept_inprogress('CSRRL05_96','greenland',[],90,[],[],[],J,'SD',0);

%%= PLOTS
% %% lets check out what some of those coefficients have been doing
% for alpha=1:round(N)
%     figure(alpha)
%     plot(thedates,slepcoffs(:,alpha))
%     datetick
%     title(sprintf('$\\alpha=%i$ Eigenvalue $=%0.3f$',alpha,V(alpha)),'interpreter','latex')
% end

%% Figure 2
% weight the first N eigentapers by G and sum
% Gfalpha(thedate,lm,alpha)
Gfalpha = zeros(length(thedates),addmoff(L),length(V)); % create matrix
for i = 1:length(thedates)
    Gfalpha(i,:,:) = G.*slepcoffs(i,:); % fill matrix
end

Gf_sum = sum(Gfalpha,3);
%% 
cax=[-1 2.5]*10^6;

figure(2)
subplot(1,2,1)
plotslep(Gf_sum(1,:)',1,[],.1,11)
caxis(cax)
titleline1 = 'Surface Density: global spherical harmonic expansion';
titleline2 = '$\sum_{l=0}^{L=90}$';
title(sprintf('%s\n%s',titleline1,titleline2),'interpreter','latex')

subplot(1,2,2)
plotplm(squeeze(dcoffs(1,:,:)),[],[],11,.1)
caxis(cax)
titleline1 = 'Surface Density: localized slepian expansion';
titleline2 = '$\\alpha$';
title(sprintf('%s\n%s',titleline1,titleline2),'interpreter','latex')

%% PLOT OF THE DIFFERENT EXPANSIONS: APRIL 2002
%some setup
% some sweet titles
titletext = {'Global Spherical Harmonics',...
    'Localized Slepian Functions',...
    'Difference'};
titlepos = [.13 0.8; .5 0.8; .87 0.8];
% some sweet subtitles
subtitletext = {'$$S(\theta,\phi)=\sum_{l=0}^{L=90}\sum_{m=-l}^{l}f_{lm}Y_{lm}(\theta,\phi)$$',...
    '$$\hat{S}(\theta,\phi)=\sum_{\alpha=1}^{\alpha=N}f_{\alpha}g_{\alpha}(\theta,\phi) $$',...
    '$$\Delta S(\theta,\phi)=S-\hat{S}$$'};
subtitlepos = [.13 0.2; .5 0.2; .87 0.2];
cax=[-1 2.5]*10^6; %coloraxis

%get positions for two subplots
figure(10)
for i = 1:3
    subplot(1,3,i)
    h = gca;
    p(i,:) = h.Position;
end
close(gcf)

%now we make the plot
f=figure(2);
bigax = axes('Parent',f,'Position',[0.1300 0.1100 0.7750 0.8150]);
axis off
for i = 1:3
    %make the axes and hide them
    ax{i} = axes('Parent',f,'Position',p(i,:));
    axis off
    %plot the data
    if i==1,plotplm(squeeze(dcoffs(1,:,:)),[],[],11,.1);
    elseif i==2,slep_lmcosi=grabanalpha(Gf_sum(1,:)',[],L,1,5);
    else 
        % difference between SH expansion and SLEP expansion
        delta_lmcosi=difflmcosi(squeeze(dcoffs(1,:,:)),slep_lmcosi);
        plotplm(delta_lmcosi,[],[],11,.1);
    end
    % set the coloraxis
    caxis(cax)
    
    set(f, 'currentaxes', bigax); 
    % some sweet titles
    text(titlepos(i,1),titlepos(i,2),titletext{i},'interpreter','latex','horizontalalignment','center')
    % some sweet subtitles
    text(subtitlepos(i,1),subtitlepos(i,2),subtitletext{i},'interpreter','latex','horizontalalignment','center')
    % set the coloraxis
    caxis(cax)
    colorbar('location','manual','position',[0.93 0.35 0.0253 0.35])
end