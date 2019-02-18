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
addpath('functions')
addpath(fullfile(dir,'functions'))
setworkspace('/Users/benjamingetraer/Documents/IndependentWork/SH_Workspace');


datenum(2012,07,07:11);

[Zdata,thetimedata,thespacelim,thelevdata] = ...
    getMerra2VAR(datenum(2012,07,07:11), 'H');

[udata] = ...
    getMerra2VAR(datenum(2012,07,07:11), 'U');
[vdata] = ...
    getMerra2VAR(datenum(2012,07,07:11), 'V');
%% Project Greenland into the image basis
[ Flon2x,Flat2y,gx,gy,bx,by,contx,conty ] = projectMERRA( Zdata, thespacelim );
azores = [Flon2x( -27.853785), Flat2y(38.471189)];
reykjavik = [Flon2x(-21.933333), Flat2y(64.133333)];

%% Eastern wind component over time
figure(1)
clf
% above 250 hPa
lvl = 22; % 17 = 500
z_slice= Zdata(:,:,lvl,1);
cax=[prctile(z_slice(:),3) prctile(z_slice(:),98)];

for t=10 %1:length(thetimedata)
    clf
    %calculation
    z_slice = imrotate(Zdata(:,:,lvl,t),90);
    js = 130/60/60*1000;
    windspeed = imrotate(sqrt(udata(:,:,lvl,t).^2 + vdata(:,:,lvl,t).^2),90);
    
    
    % plotting
    
    
    
    imagesc(z_slice,'AlphaData',~isnan(z_slice));
    hold on

    
    %     plot(gx,gy,'w-')
    plot(bx,by,'w:')
    plot(contx,conty,'k-','linewidth',1)
    
    %put a title on
    title(sprintf('Height [10^4m] of %0.0f hPa Isobar',...
        thelevdata(lvl)))
    % color
    colormap(jet)
    caxis(cax)
    %     colorbar
    %axis
    latlonaxis( xlim, ylim, thespacelim, 5 );
    % contour
    contour(z_slice./1E4,[1.03, 1.04, 1.05, 1.07, 1.08, 1.09],'k-','ShowText','on')
    
    % label high and low
    text([47 126],[71 39],'LO','fontsize',20,'horizontalalignment','center')
    text(90,68,'HI','fontsize',20,'horizontalalignment','center')
     
    handle = flowArrow(udata(:,:,lvl,t),vdata(:,:,lvl,t),0.2,2,'k');
    
    % label azores and reykjavik
    text(azores(1),azores(2),'Azores','horizontalalignment','center',...
        'BackgroundColor','w','interpreter','latex');
    text(reykjavik(1),reykjavik(2),'Reykjav\''ik','horizontalalignment','center',...
        'BackgroundColor','w','interpreter','latex');
   
    % set(gca,'xlim',[5.7531   87.8223],'ylim',[16.8935  115.3161])
    
    %legend
    legend(handle,sprintf('Windspeed: range %0.1f to %0.1f m/s',min(windspeed(:)),max(windspeed(:))));
    
     % put the date on
    text(Flon2x(-89),Flat2y(89),thetimedata(t),'color','w');
    
    % NAO data
%     NAO = dlmread(fullfile(datadir,'NAO_monthly.txt'));
%     NAOdates = datenum(NAO(:,1),NAO(:,2),repmat(15,size(NAO(:,1))));
%     NAOdata = NAO(:,3);
    
%     thisNAO = NAOdata(year(NAOdates)==2012 & month(NAOdates)==7);
    thisNAO = readDailyNAO(round(datenum(thetimedata(t),'dd-mmm-YYYY hhMM')),1);
    text(Flon2x(-89),Flat2y(87),sprintf('NAO = %0.3f',thisNAO),'color','w');
    % animation
    pause(0.5)
end
%%
lvl = 22;


figure(3)
subplot(1,3,1)

imagesc(udata(:,:,lvl,t));
colorbar
subplot(1,3,2)

imagesc(vdata(:,:,lvl,t));
colorbar
subplot(1,3,3)

js = 130/60/60*1000;
windspeed = sqrt(udata(:,:,lvl,t).^2 + vdata(:,:,lvl,t).^2);
windspeed(windspeed<=js)=0;

imagesc(imrotate(windspeed,90));
colorbar