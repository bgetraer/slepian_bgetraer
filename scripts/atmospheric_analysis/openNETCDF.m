%**************************************************************************
% Example code for opening and using netCDF files
% Merra2 Data from:
% https://disc.gsfc.nasa.gov/datasets/M2I6NPANA_V5.12.4/summary?keywords=%22MERRA-2%22
%
%**************************************************************************

% example netCDF file directory
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
merra_dir = fullfile(dir,'datafiles/MERRA2');
ex_ncdf = fullfile(merra_dir,'MERRA2_400.inst6_3d_ana_Np.20120711.nc4');

addpath('functions')

% get the file's meta-data
finfo = ncinfo(ex_ncdf);

% display the file's variables
disp('Variable List');
disp({finfo.Variables.Name}');

disp('Global Attributes');
disp({finfo.Attributes.Name}');

% Pull the data resolution attribute
disp('Spatial Resolution');
res_info = ncreadatt(ex_ncdf,'/','DataResolution'); 
disp(res_info);

% Pull the temporal range attribute
disp('Temporal Range');
timerange = ncdf_timerange(ex_ncdf);
timeRange = horzcat(timerange{1},  '  ---  ', timerange{2})
%% Pick out data of specific variables
lon_data = ncread(ex_ncdf,'lon');
lat_data = ncread(ex_ncdf,'lat');
lev_data = ncread(ex_ncdf,'lev');
time_data = ncread(ex_ncdf,'time');


temp_data = ncread(ex_ncdf,'T');
hum_data = ncread(ex_ncdf,'QV');
east_data = ncread(ex_ncdf,'U');
north_data = ncread(ex_ncdf,'V');
vector_wind_data = sqrt(east_data.^2+north_data.^2);

%% Temperature over time
figure(1)
clf
t_slice=temp_data(:,:,1,1);
cax=[prctile(t_slice(:),5) prctile(t_slice(:),95)];
for t=1:4
    t_slice = temp_data(:,:,1,t);

    imagesc(imrotate(t_slice,90))
    caxis(cax)
    colorbar
    pause(0.3)
end

%% Temperature over pressure level
close all
figure(1)
clf
t_slice=temp_data(:,:,10,1);
cax=[prctile(t_slice(:),2) prctile(t_slice(:),98)];
for p=1:42
    t_slice = temp_data(:,:,p,1);
    imagesc(imrotate(t_slice,90))
    title(sprintf('Air Temperature (K) at pressure level: %0.2f kPa',lev_data(p)/10))
    text(0,350,timerange(1),'color','w')
    caxis(cax)
    colorbar
    pause(0.1)
end
%% Specific Humidity over time (non-integrated)
figure(1)
clf
h_slice= hum_data(:,:,10,1);
cax=[prctile(h_slice(:),3) prctile(h_slice(:),98)];

for t=1:4
    %calculation
    h_slice = hum_data(:,:,10,t);
    thedate = datestr(datenum(year(timerange(1)),month(timerange(1)),...
        day(timerange(1)),hour(timerange(1)),double(time_data(t)),0),...
        'dd-mmm-yyyy HHMM');
    % plotting
    imagesc(imrotate(h_slice,90));
    text(0,350,thedate,'color','w');
    title(sprintf('Specific humidity [kg kg^{-1}] at %0.0f kPa pressure level',...
        lev_data(10)/10))
    line([144.5 288.5],[90.75,90.75],'color','w','linewidth',2)
    line([144.5 288.5],[0,0],'color','w','linewidth',2)
    line([144.5 144.5],[0,90.75],'color','w','linewidth',2)
    line([288.5 288.5],[0,90.75],'color','w','linewidth',2)
    % color
    caxis(cax)
    colorbar
    % animation
    pause(0.3)
end
%% Specific Humidity over time (integrated)
figure(1)
clf
% integrate over 1000 - 200 hPa
PL_index = lev_data<=1000 & lev_data>=200;
dP = diff(lev_data(PL_index));

hum_data(:,:,PL_index,1).*vector_wind_data(:,:,PL_index,1);
%%

for i = 1:length(dP)
    hum_data(:,:,i,1).*vector_wind_data(:,:,i,1).*dP(i)
end

h_slice=nansum(hum_data(:,:,:,1),3);
cax=[prctile(h_slice(:),5) prctile(h_slice(:),95)];

imagesc(imrotate(h_slice,90));
colorbar
%%
for t=1:4
    %calculation
    h_slice = nansum(hum_data(:,:,:,t),3);
    thedate = datestr(datenum(year(timerange(1)),month(timerange(1)),...
        day(timerange(1)),hour(timerange(1)),double(time_data(t)),0),...
        'dd-mmm-yyyy HHMM');
    % plotting
    imagesc(imrotate(h_slice,90));
    text(0,350,thedate,'color','w');
    title('Vertically integrated specific humidity [kg kg^{-1}]')
    line([144.5 288.5],[90.75,90.75],'color','w','linewidth',2)
    line([144.5 288.5],[0,0],'color','w','linewidth',2)
    line([144.5 144.5],[0,90.75],'color','w','linewidth',2)
    line([288.5 288.5],[0,90.75],'color','w','linewidth',2)
    % color
    caxis(cax)
    colorbar
    % animation
    pause(0.3)
end
%% Lat Lon gridding
figure(1)
clf
h_slice=nansum(hum_data(:,:,:,1),3);
imagesc(imrotate(h_slice,90));
text(0,350,timerange(1),'color','w');
title('Vertically integrated specific humidity [kg kg^{-1}]')
xlim = get(gca,'XLim');
xtix = linspace(xlim(1),xlim(2),9);
xtiklab = linspace(-180,180,9);
ylim = get(gca,'YLim');
ytix = linspace(ylim(1),ylim(2),7);
ytiklab = linspace(90,-90,7);
set(gca,'xtick',xtix,'xticklabel',xtiklab,...
    'ytick',ytix,'yticklabel',ytiklab,...
    'xgrid','on','ygrid','on','xminorgrid','on','yminorgrid','on',...
    'minorgridcolor','w','gridcolor','w')

%box

Fx = griddedInterpolant(xtiklab,xtix);
Fy = griddedInterpolant(flip(ytiklab),flip(ytix));

a = -90;
b = 0;
c = 30;
d = 90;

line([Fx(a) Fx(b)],[Fy(c),Fy(c)],'color','w','linewidth',2)
line([Fx(a) Fx(b)],[Fy(d),Fy(d)],'color','w','linewidth',2)
line([Fx(a) Fx(a)],[Fy(d),Fy(c)],'color','w','linewidth',2)
line([Fx(b) Fx(b)],[Fy(d),Fy(c)],'color','w','linewidth',2)
% color
colorbar

%%  wind
xwind = east_data(:,:,10,1);
ywind = north_data(:,:,10,1);
vecwind = vector_wind_data(:,:,10,1);

figure(1)
clf
subplot(1,3,1)
imagesc(imrotate(xwind,90));
title('x wind component')
subplot(1,3,2)
imagesc(imrotate(ywind,90));
title('y wind component')
subplot(1,3,3)
imagesc(imrotate(vecwind,90));
title('wind magnitude')
