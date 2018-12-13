%**************************************************************************
% Example code for opening and using netCDF files
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

[Hdata,thetimedata,thespacelim,thelevdata] = ...
    getMerra2VAR(datenum(2012,07,07:11), 'H');

%% Project Greenland into the image basis
[ Flon2x,Flat2y,gx,gy,bx,by,contx,conty ] = projectMERRA( Hdata, thespacelim );
%% Specific Humidity over time (non-integrated)
figure(1)
clf
lvl = 17;
h_slice= Hdata(:,:,lvl,1);
cax=[prctile(h_slice(:),3) prctile(h_slice(:),98)];

for t=1:20
    %calculation
    h_slice = imrotate(Hdata(:,:,lvl,t),90);
    % plotting
    imagesc(h_slice,'AlphaData',~isnan(h_slice));
    hold on
%     plot(gx,gy,'w-')
    plot(bx,by,'w:')
    plot(contx,conty,'k-','linewidth',1)
    text(Flon2x(-89),Flat2y(89),thetimedata(t),'color','w');
    title(sprintf('Specific humidity [kg kg^{-1}] at %0.0f hPa pressure level',...
        thelevdata(lvl)))
    % color
    colormap(jet)
    caxis(cax)
    colorbar
    %axis
    latlonaxis( xlim, ylim, thespacelim, 5 );
    % animation
    pause(0.3)
end
