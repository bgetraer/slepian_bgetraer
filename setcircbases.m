function [ centers,txtlocation,ID,rotation ] = setcircbases(type,r,p,rot)
%SETCIRCBASES Exports centers of rotated circular caps centered on Oceans
%or Continents. Optionally plots and labels the circles on a map.
%
%INPUT
%   type  geographic centers
%           1   major oceans (N Atl, N Pac, S Atl, S Pac, Indian)
%           2   continents (Afr, Aust, S Amer, N Amer, Asia, Eur, Grnlnd)
%               default = 1
%   r     Angular extent of the spherical cap (degrees)
%               default = 15
%   p     0     no plot
%         1     DO plot
%               default(0)
%   rot     0   do not rotate mapping longitude
%           1   rotate longitude to avoid cutting off lines
%               default(1)
%OUTPUT
%   centers     array of geographic centers in [lon;lat]
%
%   See also:
%       CAPLOC
%
% Last modified by bgetraer@princeton.edu, 11/24/2017

% set default values
defval('type',1)
defval('r',15)
defval('p',0)
defval('rot',1)

% call data from subfunction CENTERSDATA
[centers,txtlocation,ID,rotation] = centersdata(type,r);

if p==1
    % plot the bases
    if rot==1
        plotbases(r,centers,txtlocation,ID,rotation)
    elseif rot==0
        plotbases(r,centers,txtlocation,ID)
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [centers,txtlocation,ID,rotation] = centersdata(type,r)
%OCEANCENTERS Returns data related to plotting circular caps rotated over
%the major oceans or continents
%
%INPUT
%   type    1   major ocean centers
%           2   major continent center
%OUTPUT
%   vectors and matrices for locating and plotting the centers
switch type
    case 1
        % text names for labeling
        ID = {'North Atlantic' 'North Pacific' 'South Atlantic' ...
            'South Pacific' 'Indian Ocean'};
        % corresponding geographic centers in [lon;lat]
        centers = [318  190  347  225  80;...
            27  27  -27  -27  -23];
        % location for plotting labels:
        %       place above for NH, below for SH
        txtlocation = [centers(1,:); ...
            (centers(2,:)-3/2*r).*double(centers(2,:)<0) + ...
            (centers(2,:)+3/2*r).*double(centers(2,:)>0)];
        % longitudinal rotation to avoid cutoff
        rotation = 20;
    case 2
        % text names for labeling
        ID = {'Africa' 'Australia' 'S. America' 'N. America' 'Asia' ...
            'Europe' 'Greenland'};
        % corresponding geographic centers in [lon;lat]
        centers = [25  134  304  262  90  8  316;...
            8  -25  -13   43   43  45 72];
        % location for plotting labels:
        %       place above for NH, below for SH
        txtlocation = [centers(1,:); ...
            (centers(2,:)-3/2*r).*double(centers(2,:)<0) + ...
            (centers(2,:)+3/2*r).*double(centers(2,:)>0)];
        % longitudinal rotation to avoid cutoff
        rotation = 65;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotbases(r,centers,txtlocation,ID,rotation)
%PLOTBASES Create a labeled plot of circular bases centered over given
%points, possibly rotated by longitude to avoid cutting off a basis
%
%
%INPUT
%   r               radius of the bases
%   centers         matrix of centers [lon;lat]
%   txtlocation     matrix of text label locations [lon;lat]
%   ID              array of text labels
%   rotation        degrees longitude (default=360 no rotation)

defval('rotation',360)

% plot of circular bases in 'rotated' lat/lon cartesian space
xlimits = [-360+rotation rotation]; % new xlimits

% plot the continent outlines twice for 'rotation'
plotcont([0 90],[rotation -90],0,0);plotcont([rotation 90],[360 -90],0,-360);
hold on;
% plot each circular basis
for i = 1:length(centers)
    % calculate points to plot
    [circLON, circLAT] = caploc([centers(1,i) centers(2,i)],r,100,1);
    % plot and label with 'rotation'
    plot(circLON(circLON>rotation)-360, circLAT(circLON>rotation),'k.','markersize',0.5);
    plot(circLON(circLON<rotation), circLAT(circLON<rotation),'k.','markersize',0.5);
    if centers(1,i)<rotation
        text(txtlocation(1,i), txtlocation(2,i),ID{i},...
            'HorizontalAlignment','center');
    else
        text(txtlocation(1,i)-360, txtlocation(2,i),ID{i},...
            'HorizontalAlignment','center');
    end
end
axis([xlimits -90 90])
xlabel('\circ Longitude');ylabel('\circ Latitude');
end