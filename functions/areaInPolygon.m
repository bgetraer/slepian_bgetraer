function [ polygonpassindex ] = ...
    areaInPolygon( C, S, wavename, level, xp, yp, polygonX, polygonY, waveThreshindex )
%AREAINPOLYGON Creates a binary pass index for a wavelet transformation
% determined by the percentage of the wavelet support within a given
% polygon. If given, it ignores wavelets coefficients which have been
% thresholded by a second variable.
%
%   NOTE: the function runs through the wavelets by level, and does not
%   distinguish between diagonal, horizontal, or vertical detail. Thus if
%   any of the 3 detail coeficients for a level are "on," it looks at all
%   3 of them. Additionally, the area percentage threshold is hardcoded as
%   follows:
%       area threshold a wavelet of a given level must occupy within the
%       polygon (by percent of wavelet support) in order to pass:
%           highest resolution, 90%
%           top two lowest resolutions, take them all
%
% SEE ALSO:
%   WAVELEVELINDEX, INPOLYGON
%
% Last modified by bgetraer@princeton.edu, 3/12/2019

defval('waveThreshindex',ones(1,length(C)));
% empty coefficient array
C0 = zeros(1,length(C));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FIND THE WAVELETS WHICH PASS THE AREA THRESHOLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ I, ~, coefinlevel] = wavelevelINDEX( C,S );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the area threshold a wavelet of a given level must occupy within the
% polygon (by percent of wavelet support) in order to pass.
%   highest resolution, 90%
%   top two lowest resolution, take them all
areathresh = [linspace(0.9,0,level-1) 0];

% initialize pass array 
if lengthyProcessFlag('wavelet in polygon index')
    polygonpassindex = zeros(size(C));
    for i = 1:level
        fprintf('level %d \n',i)
        % if we are keeping them all anyways!
        if areathresh(i)==0
            polygonpassindex(I{i})=1;
        else
            % only look at one of the wavelet sequences
            index = I{i}(1:coefinlevel(i));
            % passes per level
            passlevel = zeros(size(index));
            for in = index
%                 fprintf('wavelet %d of %d \n',find(index==in),coefinlevel(i))
                Cmod = C0;
                % don't do the work if we don't want it anyways!
                if any([waveThreshindex(in), waveThreshindex(in+length(passlevel)),...
                       waveThreshindex(in+2*length(passlevel))])
                   % turn on the wavelet we are checking
                    Cmod(in) = 1;
                    imj = waverec2(Cmod,S,wavename);
                    % what points are in the wavelet and polygon
                    xw = xp(imj(:)~=0);
                    yw = yp(imj(:)~=0);
                    pw = inpolygon(xw(:),yw(:),polygonX,polygonY);
                    aw = sum(imj(:)~=0);
                    ainp = sum(pw(:)~=0);
                    passlevel(index==in) = (1-(aw-ainp)/aw > areathresh(i));
                else 
                    % if the wavelet is thrown out by a second index
                    passlevel(index==in) = 0;
                end
            end
            polygonpassindex(I{i}) = [passlevel, passlevel, passlevel];
        end
    end
end
polygonpassindex = polygonpassindex.*waveThreshindex;

end

