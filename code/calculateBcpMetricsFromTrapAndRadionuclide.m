function [compilation] = calculateBcpMetricsFromTrapAndRadionuclide(nLocs,...
    latsLocal,lonsLocal)

% CALCULATEBCPMETRICSFROMTRAPANDRADIONUCLIDE Calculates three metrics of the 
% the BCP (Teff 100 to 1000 m, b and z*) using the sediment trap and 
% radionuclide POC flux observational data compiled for this study.
%
%   INPUT:
%       nLocs     - number of locations
%       latsLocal - list of latitudes for interpolation
%       lonsLocal - list of longitudes for interpolation
%
%   OUTPUT:
%       compilation - local annual average values and standard
%                     deviations of b, z* and Teff 100 to 1000 m
%
%   This script uses these external functions: 
%       propagateErrorWithMCforPEeffAndTeff.m    - custom function
%       propagateErrorWithMCforMartinbAndZstar.m - custom function
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 9 April 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Definitions

% Filename declarations 
filenameMonthlyPocFlux = 'pocflux_compilation_monthly.mat';

% Define parameters
% NUM_NPP_ALGOS = 3; % CAFE, VGPM, CbPM (for calculation of PEeff)
MAX_NUM_DEPTHS_PER_PROFILE = 100;

% Load the trap and radionuclide compilation (mmol C m-2 d-1)
load(fullfile('.','data','processed',filenameMonthlyPocFlux),...
    'obsMonthlyDhAvg','obsMonthlyDhErrTot','obsMonthlyProfileDepths',...
    'obsMonthlyProfileAvg','obsMonthlyProfileErrTot')

%% Calculate PEeff and Teff

% Initialise array to offload POC flux data
% 1st dimension: 1=at zeu, 2=at zmeso
% 4th dimension: 1=avg, 2=net error 
arrayFlux = NaN(2,12,nLocs,2); 

% Extract POC flux data, with dimensions nLocs x 12
xzeu = squeeze(obsMonthlyDhAvg(:,1,:)); % mean
xzmeso = squeeze(obsMonthlyDhAvg(:,2,:)); 
xezeu = squeeze(obsMonthlyDhErrTot(:,1,:)); % total error
xezmeso = squeeze(obsMonthlyDhErrTot(:,2,:)); 

% Reshape the extracted POC flux data to dimensions 12 x nLocs, and 
% create arrayFlux to be passed to the function computing PEeff and Teff 
% from depths of 100 to 1000 m
arrayFlux(1,:,:,1) = permute(xzeu,[2,1]); 
arrayFlux(1,:,:,2) = permute(xezeu,[2,1]); 
arrayFlux(2,:,:,1) = permute(xzmeso,[2,1]);
arrayFlux(2,:,:,2) = permute(xezmeso,[2,1]);

fprintf('\nInitiate calculation of Teff...\n')
        
[teffAnnual,~] = propagateErrorWithMCforPEeffAndTeff(...
    nLocs,latsLocal,lonsLocal,[],[],arrayFlux);

fprintf('\n...done.\n')

compilation.teff100to1000.ave      = teffAnnual(:,1);
compilation.teff100to1000.stdevupp = teffAnnual(:,2);
compilation.teff100to1000.stdevlow = teffAnnual(:,3);
compilation.teff100to1000.max      = teffAnnual(:,4);
compilation.teff100to1000.min      = teffAnnual(:,5);

%% Calculate Martin's b and z*

% Initialise arrays to offload POC flux data
arrayFlux   = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,nLocs,2); % 4th dimension: 1=median, 2=net error 
arrayDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,nLocs);
for iLoc = 1:nLocs
    for iMonth = 1:12
        izEuAndBelow = obsMonthlyProfileDepths(:,iMonth,iLoc) >= 100 &...
            ~isnan(obsMonthlyProfileDepths(:,iMonth,iLoc));
        nEntries = sum(izEuAndBelow == 1);
        arrayFlux(1:nEntries,iMonth,iLoc,1) = obsMonthlyProfileAvg(izEuAndBelow,iMonth,iLoc);
        arrayFlux(1:nEntries,iMonth,iLoc,2) = obsMonthlyProfileErrTot(izEuAndBelow,iMonth,iLoc);
        arrayDepths(1:nEntries,iMonth,iLoc) = obsMonthlyProfileDepths(izEuAndBelow,iMonth,iLoc);
    end
end

fprintf('\nInitiate calculation of b and z*...\n')

[martinbAnnual,zstarAnnual,martinbMonthly_gof,zstarMonthly_gof] =... 
    propagateErrorWithMCforMartinbAndZstar(nLocs,arrayDepths,arrayFlux);

fprintf('\n...done.\n')

compilation.martinb.ave      = martinbAnnual(:,1);
compilation.martinb.stdevupp = martinbAnnual(:,2);
compilation.martinb.stdevlow = martinbAnnual(:,3);
compilation.martinb.max      = martinbAnnual(:,4);
compilation.martinb.min      = martinbAnnual(:,5);
compilation.martinb.gof      = martinbMonthly_gof;

compilation.zstar.ave        = zstarAnnual(:,1);
compilation.zstar.stdevupp   = zstarAnnual(:,2);
compilation.zstar.stdevlow   = zstarAnnual(:,3);
compilation.zstar.max        = zstarAnnual(:,4);
compilation.zstar.min        = zstarAnnual(:,5);
compilation.zstar.gof        = zstarMonthly_gof;

%% Save metrics

save(fullfile('.','data','processed','bcpmetrics_compilation.mat'),'compilation*')

end
