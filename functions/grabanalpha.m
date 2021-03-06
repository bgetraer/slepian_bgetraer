function [ lmcosi_mat,lmcosi_sum,alphaindex ] = grabanalpha(G,N,L,alphaindex,plot)
%GRABANALPHA
%   "Grabs" the lm coefficients from a given order alpha of a Slepian basis
%   or bases given by (lm)X(alpha)X(basis) matrix G, such as the
%   output of: GLMALPHA GLMALPHAPTO or SLEPCIRCBASES
%
%INPUTS
%   G           The (lm)X(alpha)X(basis) square matrix of coefficients
%                   or 'greenland'
%   N           Shannon numbers of the alphas
%   L           Bandwidth of the alphas
%   alphaindex  Index of alphas of length size(G,3), or a single alpha to
%                   be grabbed for all bases. Default = 
%   plot        Plot the output?
%                   0   NO (default)
%                   1   Mollweide projection
%                   2   3-D sphere, no topography
%                   3   Flat rectangular projection for entire globe
%                   4, or ~=0 and 'greenland' was asked for
%                           plot pretty for greenland
%                   5   like 4 but default colormap
%
%OUTPUTS
%   lmcosi_mat      A stack of standard spherical harmonic coefficients
%                   matrices corresponding to the "grabbed" alphas.
%   lmcosi_sum      A single standard spherical harmonic coefficients
%                   matrix corresponding to the SUM of the "grabbed"
%                   alphas.
%   alphaindex      Returns the vector of alpha values, ie the rank of the Slepian
%                   function, helpful when alpha is grabbed randomly.
%
%       See also:
%           SLEPCIRCBASES GLMALPHAPTO GLMALPHA ADDMON PLM2XYZ PLOTPLM
%
%Last modified by bgetraer@princeton.edu, 12/04/2017

%did we ask for the G of Greenland?
if ischar(G) && unique(G=='greenland')==1
    defval('L',90)
    %is this G already there?
    if exist(fullfile(getenv('IFILES'),'KERNELC',sprintf('WREG-greenland-%i.mat',L)),'file')~=2
        error('KernelC File for requested bandlimit has not been created yet.');
    end
    %it is? load it up (could take a bit)
    G = glmalpha('greenland',L);
    %did we want to plot? well lets do it pretty for greenland
    if plot~=0,plot=4;end
end

%if no alphaindex is given, randomize
defval('N',8)
defval('alphaindex',max(randi(round(N*1/2),1,size(G,3)),1));
defval('plot',0)

%if one alphaindex is given, return only that
if length(alphaindex)==1 && length(alphaindex)<size(G,3)
    alphaindex = repmat(alphaindex,1,size(G,3));
end

% Create blank LMCOSI matrix
[~,~,~,blank,~,~,~,~,~,ronm]=addmon(L);
lmcosi_mat = zeros([size(blank) size(G,3)]) + blank;

for j = 1:size(G,3)
    % Create the coefficient blanks
    cosi=blank(:,3:4);
    % grab the coefficients of an alpha eigentaper and
    % re-index them in lmcosi format
    cosi(ronm)=G(:,alphaindex(j),j);
    % Add them to the full matrix
    lmcosi_mat(:,3:4,j)=cosi;
end
% Create 2d matrix of summed coefficients for all grabbed alphas.
lmcosi_sum = [blank(:,1:2) sum(lmcosi_mat(:,3:4,:),3)];

% Did we want a plot?
switch plot
    case 0
    case 1
        plotplm(lmcosi_sum,[],[],7,1)
    case 2
        plotplm(lmcosi_sum,[],[],2,.1)
    case 3
        plotplm(lmcosi_sum,[],[],9,.1)
    case 4
        plotplm(lmcosi_sum,[],[],10,.1)
    case 5
        plotplm(lmcosi_sum,[],[],11,.2)
    otherwise
        error('Not a valid plotting method')
end
%title?
% if length(unique(alphaindex))==1 && plot~=0
%     title(strcat('$\alpha=$',sprintf('%i',unique(alphaindex))),'interpreter','latex');
% end
end