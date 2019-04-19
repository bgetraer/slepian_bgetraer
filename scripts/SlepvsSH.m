%%*************************************************************************
% Comparison of Spherical Harmonic Eigenfunctions to Slepian Eigentapers in
% a circular basis and in the Greenland Geographic bases.
% Generates an illustrative figure.
% bgetraer
%% SETUP WORKSPACE
setworkspace('/Users/benjamingetraer/Documents/IndependentWork/SH_Workspace');

% demonstrate SH functions for a few orders
lindex = [1,1,2];
mindex = [0,1,1];

%choose a circular region (ie N Atl)
r = 15;
L = 20;
[Gcirc,N] = slepcircbases(r,L,1);
Gcirc = Gcirc(:,:,1);
%find eigentapers corresponding to SH eigenfunctions
%imtap corresponds to order l, mtap corresponds to degree +-m
%for axisymmetric slepian bases
[~,~,~,~,~,~,mtap,imtap]=glmalphapto([],L,[],[],[]);
for j = 1:length(lindex)
    alphaindex(1:2,j) = find(imtap==lindex(j)&(mtap==mindex(j)|mtap==-mindex(j)))
end
%the Greenland L=90 corresponding eigentapers
greenlandalpha = [1 2 7 5];

%% NOW PLOT IT ALL
% choose 3d plotting projection
plottype=2;
% and colormap
cmap = 'bluewhitered([],1)';


%create dummy figure and pull locations from there
figure(100);
for i = 1:9
    subplot(3,3,i)
    h = gca;
    p1(i,:) = h.Position;
end
close(gcf)

p2 = p1([5,6,8],:)+[0.065 0 0 0]; % axes for the second double figures
p1([5,6,8],:) = p1([5,6,8],:)-[0.065 0 0 0]; % axes for the first double figures
p = [p1;p2]; % axes for all figures.

% generate blank figure and axes
f = figure(2)
clf
for j = 1:size(p,1)
    ax{j} = axes('Parent',f,'Position',p(j,:));
    axis off
end

% fill the axes
for i = 1:3
    %first row
    set(f, 'currentaxes', ax{i}); %axes 1 2 and 3
    shpower(lindex(i),mindex(i),plottype)
    if plottype==2,view(48,18);end

    %second row
    set(f, 'currentaxes', ax{i+3}); %axes 4 5 and 6
    grabanalpha(Gcirc,N,L,alphaindex(1,i),plottype);
    if plottype==2,view(48,18);end
    
    %second row, double figures
    if i<3
        set(f, 'currentaxes', ax{i+9}); %axes 10 and 11
        grabanalpha(Gcirc,N,L,alphaindex(2,i+1),plottype);
        if plottype==2,view(48,18);end
    end
    
    %third row
    set(f, 'currentaxes', ax{i+6}); %axes 7 8 and 9
    grabanalpha('greenland',[],[],greenlandalpha(i),plottype);
    
    %third row, double figures
    if i==3
        set(f, 'currentaxes', ax{i+9}); %axis 12
        grabanalpha('greenland',[],[],greenlandalpha(i+1),plottype);
    end
end

% colorbar
ax_big = axes('Parent',f,'Position',[0.1300 0.1100 0.7750 0.8150]);
axis off
caxis([-1,1]);
cb = colorbar('location','manual');
cb.Position = [0.91 0.08 0.03 0.3];
cb.Ticks = -1:0.5:1;
colormap(eval(cmap));

% text labels
circlelabel = sprintf('Circular basis of $r=%i$, $L=%i$',r,L);
greenlandlabel = 'Geographic basis of $L=60$';

text(0.5,0.65,circlelabel,'horizontalalignment','center',...
    'interpreter','latex','fontsize',12);
text(0.5,0.32,greenlandlabel,'horizontalalignment','center',...
    'interpreter','latex','fontsize',12);
% Save figure
figurefiledir = '/Users/benjamingetraer/Documents/JuniorPaper/JP01/Figures/SLEPexamples/SlepvSH'
print('-bestfit',figurefiledir,'-dpdf','-r300')