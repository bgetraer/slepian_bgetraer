%**************************************************************************
% Difference Maps in SH, SLEP, and MODELED
% Last modified by bgetraer@princeton.edu, 1/8/2017
%**************************************************************************
addpath('/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/functions')
setworkspace('/Users/benjamingetraer/Documents/JuniorPaper/SH_Workspace');
datadir = '/Users/benjamingetraer/Documents/JuniorPaper/slepian_bgetraer/datafiles';
% Get 'thedates','ESTtotal','ESTtotalresid','total','alphavarall' from GREENLAND60.m
load(fullfile(datadir,'Greenland60data'));
N=20;%for 0.5 buffered greenland L=60

[dcoffs,~,ddates] = ...
    grace2plmt_inprogress('CSR','RL05','SD',0,60);
validrange = monthnum(1,2003,ddates):length(ddates);
dcoffs = dcoffs(validrange,:,:);
ddates = ddates(validrange);


%date domains:
Harig2013 = 1:monthnum(6,2013,thedates);
Getraer2018 = monthnum(6,2013,thedates):length(thedates);
before = 1:monthnum(6,2012,thedates);
during = monthnum(7,2012,thedates):monthnum(8,2013,thedates);
after = monthnum(9,2013,thedates):length(thedates);
all = 1:length(thedates);
without = [before after];

% colormap
cmap = 'bluewhitered([],1)';

%% SH DIFF MAP January to September 2012
y1=2012;y2=2012;
m1 = monthnum(9,y1,thedates);
m2 = monthnum(1,y2,thedates);
sh1 = squeeze(dcoffs(m1,:,:));
sh2 = squeeze(dcoffs(m2,:,:));

diff2012SH=difflmcosi(cat(3,sh1,sh2));
diff2012SLEP = (slepcoffs(m1,1:N)-slepcoffs(m2,1:N));

f =figure(1);
clf

p = [0   0.1100    0.5    0.8150];
ax = axes(f,'position',p);
axis off
hold on 
plotplm(diff2012SH,[],[],2,0.2);
% Here I improvise a grid
sphlatlon(10,10);
view(45,90)
cax=[-1500 1500];
colormap(eval(cmap))
caxis(cax);
title('Unfiltered','interpreter','latex')
set(gca,'fontsize',12)

p = [0.4    0.1100    0.5    0.8150];
ax2 = axes(f,'position',p);
axis off
ax1 = copyobj(ax2,f);
axis off

hold on 
plotplm(plmfilt(diff2012SH,30),[],[],2,0.2);
% Here I improvise a grid
sphlatlon(10,10)
view(45,90)
caxis(cax);
colormap(eval(cmap))

title('Filtered','interpreter','latex')
set(gca,'fontsize',12)


set(f, 'currentaxes', ax2);
axis off
caxis(cax);
colormap(eval(cmap))
c2 = colorbar;
c2.Label.String = 'kg per m$^2$';
c2.Position = [0.85 0.2667 0.0299 0.5019];
c2.Label.Interpreter = 'latex';
set(gca,'fontsize',12)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT COMPARISON OF THE DIFFERENT EXPANSIONS: APRIL 2002
%some setup
% some sweet titles
titletext = {'Global Spherical Harmonics',...
    'Localized Slepian Functions',...
    'Difference'};
titlepos = [.13 0.8; .5 0.8; .87 0.8];
% some sweet subtitles
subtitletext = {'$$S_{L}(\theta,\phi)=\sum_{l=0}^{L=60}\sum_{m=-l}^{l}f_{lm}Y_{lm}(\theta,\phi)$$',...
    '$$\hat{S}_{N}(\theta,\phi)=\sum_{\alpha=1}^{N=20}f_{\alpha}g_{\alpha}(\theta,\phi) $$',...
    '$$\Delta S(\theta,\phi)=S_{L}-\hat{S}_{N}$$'};
subtitlepos = [.13 0.2; .5 0.2; .87 0.2];
cax=[-750 500]; %coloraxis

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
clf
clear ax
bigax = axes('Parent',f,'Position',[0.1300 0.1100 0.7750 0.8150]);
axis off
for i = 1:3
    %make the axes and hide them
    ax{i} = axes('Parent',f,'Position',p(i,:));
    axis off
    %plot the data
    if i==1,plotplm(diff2012SH,[],[],11,.1);
    elseif i==2,diff2012SLEP_lmcosi=grabanalpha(G(:,1:N)*diff2012SLEP',[],60,1,5);
    else 
        % difference between SH expansion and SLEP expansion
        dim = min([size(diff2012SH,1) size(diff2012SLEP_lmcosi,1)]);
        delta_lmcosi=difflmcosi(cat(3,diff2012SH(1:dim,1:4),diff2012SLEP_lmcosi(1:dim,1:4)));
        plotplm(delta_lmcosi,[],[],11,.1);
    end
    % set the coloraxis
    caxis(cax)
    colormap(eval(cmap))
    
    set(f, 'currentaxes', bigax);
    % some sweet titles
    text(titlepos(i,1),titlepos(i,2),titletext{i},'interpreter','latex','horizontalalignment','center','fontsize',12)
    % some sweet subtitles
    text(subtitlepos(i,1),subtitlepos(i,2),subtitletext{i},'interpreter','latex','horizontalalignment','center','fontsize',12)
end
% set the coloraxis
caxis(cax)
colormap(eval(cmap))
c2 = colorbar;
c2.Label.String = 'kg per m$^2$';
c2.Position = [0.93 0.35 0.0253 0.35];
c2.Label.Interpreter = 'latex';
set(gca,'fontsize',12)
%% and SLEP difference maps
y1=2012;y2=2012;
m1 = monthnum(9,y1,thedates);
m2 = monthnum(1,y2,thedates);
sh1 = squeeze(dcoffs(m1,:,:));
sh2 = squeeze(dcoffs(m2,:,:));
diff_coff = (slepcoffs(m1,1:N)-slepcoffs(m2,1:N));
scaling = 1;  %kg/m^2=mmwe=1/10cmwe

figure(1)
clf
plotplm(difflmcosi(cat(3,sh1,sh2),scaling),[],[],10,0.2);
cax = [-100 150];
caxis(cax);
colormap(eval(cmap))
colorbar

figure(2)
clf
plotslep(G(:,1:N)*diff_coff'*scaling,1,[],2,5);
caxis(cax);
colormap(eval(cmap))
colorbar

%%
ESTdiff_coff = (ESTsignal(m1,1:N)-ESTsignal(m2,1:N));
figure(3)
clf
plotslep(G(:,1:N)*ESTdiff_coff'*scaling,1,[],2,5);
caxis(cax);
colormap(eval(cmap))
colorbar

%% lets check out what some of those coefficients have been doing
% for alpha=1:round(N)
%     figure(alpha)
%     plot(thedates,slepcoffs(:,alpha))
%     datetick
%     title(sprintf('$\\alpha=%i$ Eigenvalue $=%0.3f$',alpha,V(alpha)),'interpreter','latex')
% end

%% Figure 2
% weight the first N eigentapers by G and sum
% Gfalpha(thedate,lm,alpha)

Gf_sum = zeros(length(thedates),addmoff(60)); % create matrix
for i = 1:length(thedates)
    Gf_sum(i,:) = G(:,1:N)*slepcoffs(i,1:N)'; % fill matrix
end

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
titleline2 = '\\ $\alpha$';
title(sprintf('\\begin{tabular}{c} %s %s \\end{tabular}',titleline1,titleline2),'interpreter','latex')

