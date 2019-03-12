%**************************************************************************
% GREENLAND 60 - expand spherical harmonics of bandwidth L=60 into the
% bandwidth L=60 slepian basis for Greenland, save relevant files.
%
% NOTE: this script makes use of the functions GRACE2SLEPT_INPROGRESS and
% GRACE2PLMT_INPROGRESS which are modified to handle the CSR RL05 96X96
% data product.
%
% Greenland60 is a script that saves
%       1) the slepian expansion coefficients in kg/m^2
%       2) the area integrated slepian eigentapers on a unit sphere
%       3) the slepian modeled signal of the GRACE time-series scaled to
%       an Earth sized spheroid in Gt. 
%
% The new files are located in the datadir directory under Greenland60data
% with the names
%       1) slepcoffs60
%       2) thedates
%       3) G
%       4) basis60INT.m
%       5) intmass60.m
%
% Last modified by bgetraer@princeton.edu, 1/4/2017
%**************************************************************************
addpath('/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/functions')
setworkspace('/Users/benjamingetraer/Documents/IndependentWork/SH_Workspace');
datadir = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/datafiles';


%% PART 1: Coefficients to timeseries and model calculations 
% load the full expansion coefficients
L = 60;
[slepcoffs,slep_calerrors,thedates,TH,G,CC,V,N] = ...
    grace2slept('CSRRL05_60','greenland',0.5,L,[],[],[],[],'SD',0);

% prep for residual analysis 
validrange = monthnum(1,2003,thedates):length(thedates); % crop off the first 7 months of untrustworthy data: 
slepcoffs = slepcoffs(validrange,:);
thedates = thedates(validrange);
% look for periodic signals
annual = 365.2422; % solar year in days
semiannual = annual/2; % in days
% execute the function
[ESTsignal,ESTresid,ftests,extravalues,total,alphavarall,totalparams,...
    totalparamerrors,totalfit,functionintegrals,alphavar] = ...
        slept2resid(slepcoffs,thedates,[3 annual semiannual],[],[],CC,TH);

% these are the modeled trends
ESTtotal = functionintegrals*ESTsignal(:,1:round(N))';
ESTtotalresid = functionintegrals*(ESTresid(:,1:round(N)))';
% these are the coefficient variancess
thevars = alphavar./(functionintegrals.^2);
% now correct everything to their own mean
ESTtotal = ESTtotal-mean(ESTtotal);
ESTtotalresid = ESTtotalresid-mean(ESTtotalresid);
total = total-mean(total);
% put em in columns
ESTtotal = ESTtotal(:);
ESTtotalresid = ESTtotalresid(:);
total = total(:);
% and export
save(fullfile(datadir,'Greenland60data'),'ESTsignal','ESTresid',...
    'ftests','extravalues',...
    'thedates','ESTtotal','ESTtotalresid','total','alphavarall',...
    'functionintegrals','alphavar','slepcoffs','G')

%% Part 2: a basic example of calculating the total mass from the SUMMED coefficients.
% expand into the first N eigentapers of the Greenland 60x60 basis
L = 60; % bandwidth of the basis into which we want to expand
J = 'N';    % how many eigentapers we want
[slepcoffs60,~,thedates,~,G,~,~,~] = ...
    grace2slept_inprogress('CSRRL05_60','greenland',0.5,L,[],[],[],J,'SD',0);

% The important bits in going from coefficients to integrated mass
load(sprintf('%s/intslep_greenland_%i_N.mat',datadir,L)) %load basis60INT
% basis60INT=integratebasis(G,'greenland'); %integrated contribution to a unit sphere
intarea = basis60INT*4*pi*fralmanac('a_EGM96','Earth')^2*10^(-12); %scale area on a unit sphere to Earth, scale kg/m^2 to Gt/m^2

intmass60 = intarea*slepcoffs60'; %the unaltered slepian time-series

% crop off the first 7 months of untrustworthy data:
validrange = monthnum(1,2003,thedates):length(thedates);
slepcoffs60 = slepcoffs60(validrange,:);
thedates = thedates(validrange);
intmass60 = intmass60(validrange);

% center the integrated mass timeseries around its mean
intmass60 = intmass60-mean(intmass60);

save(fullfile(datadir,'Greenland60data'),'slepcoffs60','thedates','G','basis60INT','intmass60')