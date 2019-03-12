%**************************************************************************
% Code for figure illustrating the polar jet stream
% Merra2 Data from:
% https://disc.gsfc.nasa.gov/datasets/M2I6NPANA_V5.12.4/summary?keywords=%22MERRA-2%22
%
%
% bgetraer@princeton.edu
%**************************************************************************
addpath('/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/functions')
datadir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/datafiles';
addpath(datadir)

% example netCDF file directory
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
merra_dir = fullfile(dir,'datafiles/MERRA2');
ex_ncdf = fullfile(merra_dir,'inst6_3d_ana_Np_SUBSET.20120707.nc4');
addpath(fullfile(dir,'scripts/atmospheric_analysis/functions'))
addpath(fullfile(dir,'functions'))
setworkspace('/Users/benjamingetraer/Documents/IndependentWork/SH_Workspace');


d = datenum(2012,07,11);

[Zdata,thetimedata,thespacelim,thelevdata] = ...
    getMerra2VAR(d, 'H');

[tdata] = ...
    getMerra2VAR(d, 'T');

%% DETERMINE WHICH POINTS TO EVALUATE AT, WHICH ONES ARE IN GREENLAND

[glatlon] = greenland(10,0.5);
glon = glatlon(:,1);
glat = glatlon(:,2);

lonlm = minmax(glon);
latlm = minmax(glat);
-76,59,-9,85
fprintf('%i %i %i %i',min(glon),max(glon),min(glat),max(glat))
gridlon = linspace(min(glon),max(glon),50);
gridlat = linspace(min(glat),max(glat),50);
[X,Y] = meshgrid(gridlon,gridlat);

figure(2)
clf
hold on
axis([lonlm latlm]);
plot(X,Y,'k.')
plot(glatlon(:,1),glatlon(:,2),'r','linewidth',2)
valid = inpolygon(X,Y,glon,glat);
plot(X(valid),Y(valid),'b.')

%%



freeze = 273.15; % Kelvin
t_slice = floorData(tdata(:,:,:,1));
B = imgaussfilt(t_slice,1);

t_slice = B<=freeze;

figure(1)
clf
% above 250 hPa
lvl = 1; % 17 = 500
% t_slice= tdata(:,:,lvl,1);
% cax=[prctile(t_slice(:),3) prctile(t_slice(:),98)];

    %calculation
    t_slice = imrotate(t_slice,90);
    js = 130/60/60*1000;
    
    
    % plotting
    imagesc(t_slice,'AlphaData',~isnan(t_slice));
    hold on
    
    
    %     plot(gx,gy,'w-')
    plot(bx,by,'w:')
    plot(contx,conty,'k-','linewidth',1)
    
    %put a title on
    title(sprintf('Height [10^4m] of %0.0f hPa Isobar',...
        thelevdata(lvl)))
    % color
    colormap(jet)
%     caxis(cax)
        colorbar
    %axis
    latlonaxis( xlim, ylim, thespacelim, 5 );
   
    
     % put the date on
    text(Flon2x(-89),Flat2y(89),thetimedata(1),'color','w');
    
 
   