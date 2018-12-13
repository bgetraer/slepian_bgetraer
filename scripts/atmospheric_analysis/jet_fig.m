%**************************************************************************
% Code for figure illustrating the polar jet stream
% Merra2 Data from:
% https://disc.gsfc.nasa.gov/datasets/M2I6NPANA_V5.12.4/summary?keywords=%22MERRA-2%22
%
%**************************************************************************

% example netCDF file directory
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
merra_dir = fullfile(dir,'datafiles/MERRA2');
ex_ncdf = fullfile(merra_dir,'inst6_3d_ana_Np_SUBSET.20120707.nc4');
addpath('functions')
addpath(fullfile(dir,'functions'))
setworkspace('/Users/benjamingetraer/Documents/IndependentWork/SH_Workspace');


datenum(2012,07,07:11);

[Zdata,thetimedata,thespacelim,thelevdata] = ...
    getMerra2VAR(datenum(2012,07,07:11), 'U');
[Zdata] = ...
    getMerra2VAR(datenum(2012,07,07:11), 'H');
%% Project Greenland into the image basis
[ Flon2x,Flat2y,gx,gy,bx,by,contx,conty ] = projectMERRA( Zdata, thespacelim );
azores = [Flon2x( -27.853785), Flat2y(38.471189)];
reykjavik = [Flon2x(-21.933333), Flat2y(64.133333)];

%% Eastern wind component over time 
figure(1)
clf
% above 250 hPa
lvl = 22;
h_slice= Zdata(:,:,lvl,1);
cax=[prctile(h_slice(:),3) prctile(h_slice(:),98)];

for t=10
    %calculation
    h_slice = imrotate(Zdata(:,:,lvl,t),90);
    % plotting
    imagesc(h_slice,'AlphaData',~isnan(h_slice));
    hold on
%     plot(gx,gy,'w-')
%     plot(bx,by,'w:')
    plot(contx,conty,'k-','linewidth',1)
    % put the date on
    text(Flon2x(-89),Flat2y(89),thetimedata(t),'color','w');
    %put a title on
    title(sprintf('Height [m] of %0.0f hPa Isobar',...
        thelevdata(lvl))) 
    % color
    colormap(jet)
    caxis(cax)
%     colorbar
    %axis
    latlonaxis( xlim, ylim, thespacelim, 5 );
    % contour
    contour(h_slice./1E4,[1.03, 1.04, 1.05, 1.07, 1.08, 1.09],'k-','ShowText','on')
    
    % label high and low
    text([47 126],[71 39],'LO','fontsize',20,'horizontalalignment','center')
    text(90,68,'HI','fontsize',20,'horizontalalignment','center')
    text(azores(1),azores(2),'Azores')
    text(azores(1),azores(2),'Reykjavik')
%     contour(h_slice,[10550 10550],'k--','ShowText','on')   
    % animation
%     pause(0.3)
end
