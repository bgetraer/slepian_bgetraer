%% GREENLAND SUBREGIONS
dir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer';
datadir = fullfile(dir,'datafiles/');

[importdata] = dlmread(fullfile(datadir,'gr_drainage.csv'),' ',7);

subreg = [importdata(:,3) importdata(:,9) importdata(:,14)];

figure(6)
clf
hold on 
colour = {'b','r','k','g'};
nw = [1,2,18,19];
ne = [3:5];
se = [6:12];
sw = [13:17];
regions = {nw,ne,se,sw};
I = unique(subreg(:,1));

LON = cell(size(regions));
LAT = cell(size(regions));

for j = 1:length(regions)
    thisLon = [];
    thisLat = [];
    for i=1:length(regions{j})
        lat = subreg(subreg(:,1)==I(regions{j}(i)),2);
        lon = subreg(subreg(:,1)==I(regions{j}(i)),3);
        
        [thisLon, thisLat] = polybool('union',lon,lat,thisLon,thisLat);
    end
    LON{j} = thisLon;
    LAT{j} = thisLat;
end

% save(fullfile(datadir,'subregions'),'LON','LAT')

%% LOAD FILES
subregions = load(fullfile(datadir,'subregions'),'LON','LAT');
figure(1)
clf
hold on
for i = 1:4
    fill(subregions.LON{i}, subregions.LAT{i},colour{i})
end

gxy = greenland(10);
buff= greenland(10,0.5);
%% EXTEND TO GREENLAND COAST AND BUFFER
figure(2)
clf
hold on

for i = 2
%     fill(subregions.LON{i}, subregions.LAT{i},colour{i})
end


% % NW
% a = 1;
% b = 21050;
% c = 2230;
% d = 3000;
% e = 2430;
% f = 3160;
% [Sparselon,Sparselat] = reducem(subregions.LON{1}(a:b),subregions.LAT{1}(a:b),0.01);
% nw = [[Sparselon; gxy(c:d,1)],...
%     [Sparselat;gxy(c:d,2)]];
% nwb = [[Sparselon; buff(e:f,1)],...
%     [Sparselat;buff(e:f,2)]];
% 
% fill(nwb(:,1), nwb(:,2),'k')

a = 1250;
b = 26940;
c = 1;
d = 615;
cc = 3000;
dd = 3150;
e = 1;
f = 350;
ee = 3160;
ff = 3250;
[Sparselon,Sparselat] = reducem(subregions.LON{2}(a:b),subregions.LAT{2}(a:b),0.01);
ne = [[gxy([c:d],1);Sparselon; gxy([cc:dd],1)],...
    [gxy([c:d],2);Sparselat;gxy([cc:dd],2)]];
neb = [[buff([e:f],1);Sparselon; buff([ee:ff],1)],...
    [buff([e:f],2);Sparselat;buff([ee:ff],2)]];

fill(neb(:,1), neb(:,2),'r')
% fill(subregions.LON{i}(a:b), subregions.LAT{i}(a:b),colour{i})

a = 5040;
b = 39089;
aa = 1;
bb = 1300;
A = [a:b aa:bb ];
c = 615;
d = 1480;
e = 350;
f = 1380;
[Sparselon,Sparselat] = reducem(subregions.LON{3}(A),subregions.LAT{3}(A),0.01);
se = [[gxy([c:d],1);Sparselon],...
    [gxy([c:d],2);Sparselat]];
seb = [[buff([e:f],1);Sparselon],...
    [buff([e:f],2);Sparselat]];
fill(seb(:,1), seb(:,2),'y')
% fill(subregions.LON{3}(A), subregions.LAT{3}(A),colour{i+1})


% a = 1;
% b = 24800;
% aa = 26500;
% bb = 32468;
% A = [aa:bb a:b];
% c = 1480;
% d = 2230;
% e = 1380;
% f = 2430;
% [Sparselon,Sparselat,~,tol] = reducem(subregions.LON{4}(A),subregions.LAT{4}(A),0.01);
% sw = [[gxy([c:d],1);Sparselon],...
%     [gxy([c:d],2);Sparselat]];
% swb = [[buff([e:f],1);swSparselon],...
%     [buff([e:f],2);swSparselat]];
% fill(swb(:,1), swb(:,2),'g')

%%

LONGL = {nw(:,1), ne(:,1), se(:,1), sw(:,1)};
LATGL = {nw(:,2), ne(:,2), se(:,2), sw(:,2)};
LONBUF = {nwb(:,1), neb(:,1), seb(:,1), swb(:,1)};
LATBUF = {nwb(:,2), neb(:,2), seb(:,2), swb(:,2)};

figure(3)
clf; hold on
for i = 1:4
%     fill(LONBUF{i}, LATBUF{i},colour{i},'facealpha',0.5)
    fill(LONGL{i}, LATGL{i},colour{end-i+1},'facealpha',0.5)
end

%% FIND CENTER
figure(1)
clf
hold on
for i = 1:4
    fill(LON{i}, LAT{i},colour{i})
%     fill(LONGL{i}, LATGL{i},colour{end-i+1})
    plot(322.42,72.57,'*')
end
center = [322.42,72.57];
%% SAVE
save(fullfile(datadir,'subregions'),'LON','LAT','LONGL','LATGL',...
    'LONBUF','LATBUF','center');
