function [uvp] = calculateBcpMetricsFromUvp(nLocs,latsLocal,lonsLocal)

% CALCULATEBCPMETRICSFROMUVP Calculates four metrics of the BCP (PEeff, 
% Teff 100 to 1000 m, b and z*) using the UVP-5 derived POC flux data
% calculated using the model of Bisson et al. (2022).
%
%   INPUT: 
%       nLocs     - number of locations
%       latsLocal - list of latitudes for interpolation
%       lonsLocal - list of longitudes for interpolation
%
%   OUTPUT:
%       uvp - local annual average values and standard
%             deviations of b, z* and Teff 100 to 1000 m
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
filenameUvpProcessedDataset45sc = 'pocflux_bisson_45sc_monthly_and_annual.mat';

% Define parameters
NUM_NPP_ALGOS = 3; % CAFE, VGPM, CbPM

% Load the UVP5 dataset 
load(fullfile('.','data','processed','UVP5',filenameUvpProcessedDataset45sc),...
    'uvpMonthlyFlux','targetDepths')

% Get index to depths where BCP metrics will be calculated
[~,iz1000] = min(abs(targetDepths-1000));
%find(targetDepths>=1000,1); 
% iz100 = find(targetDepths>=100,1); 
[~,iz100] = min(abs(targetDepths-100));

%% Calculate PEeff and Teff

% Initialise array to offload POC flux data
% 1st dimension: 1=at zeu, 2=at zmeso
% 4th dimension: 1=avg, 2=net error 
arrayFlux = NaN(2,12,nLocs,2); 
arrayFlux(1,:,:,1) = uvpMonthlyFlux(iz100,:,:,1); % avg at 100 m
arrayFlux(1,:,:,2) = uvpMonthlyFlux(iz100,:,:,2); % net error at 100 m
arrayFlux(2,:,:,1) = uvpMonthlyFlux(iz1000,:,:,1); % avg at 1000 m
arrayFlux(2,:,:,2) = uvpMonthlyFlux(iz1000,:,:,2); % net error at 1000 m

fprintf('\nInitiate calculation of Teff...\n')

[teffAnnual] = propagateErrorWithMCforPEeffAndTeff(...
    nLocs,0,latsLocal,lonsLocal,arrayFlux);

fprintf('\n...done.\n')

uvp.teff100to1000.ave      = teffAnnual(:,1);
uvp.teff100to1000.stdevupp = teffAnnual(:,2);
uvp.teff100to1000.stdevlow = teffAnnual(:,3);
uvp.teff100to1000.max      = teffAnnual(:,4);
uvp.teff100to1000.min      = teffAnnual(:,5);

%% Calculate Martin's b and z*

% Initialise arrays to offload POC flux data
nDepths     = size(uvpMonthlyFlux(iz100:end,:,:,1),1); 
arrayFlux   = NaN(nDepths,12,nLocs,2); % 4th dimension: 1=median, 2=net error
arrayDepths = NaN(nDepths,12,nLocs);
for iLoc = 1:nLocs
    for iMonth = 1:12
        arrayFlux(:,iMonth,iLoc,1) = uvpMonthlyFlux(iz100:end,iMonth,iLoc,1);
        arrayFlux(:,iMonth,iLoc,2) = uvpMonthlyFlux(iz100:end,iMonth,iLoc,2);
        arrayDepths(:,iMonth,iLoc) = targetDepths(iz100:end);
    end
end

fprintf('\nInitiate calculation of b and z*...\n')

[martinbAnnual,zstarAnnual,martinbMonthly_gof,zstarMonthly_gof] =... 
    propagateErrorWithMCforMartinbAndZstar(nLocs,arrayDepths,arrayFlux);

fprintf('\n...done.\n')

uvp.martinb.ave      = martinbAnnual(:,1);
uvp.martinb.stdevupp = martinbAnnual(:,2);
uvp.martinb.stdevlow = martinbAnnual(:,3);
uvp.martinb.max      = martinbAnnual(:,4);
uvp.martinb.min      = martinbAnnual(:,5);
uvp.martinb.gof      = martinbMonthly_gof;

uvp.zstar.ave        = zstarAnnual(:,1);
uvp.zstar.stdevupp   = zstarAnnual(:,2);
uvp.zstar.stdevlow   = zstarAnnual(:,3);
uvp.zstar.max        = zstarAnnual(:,4);
uvp.zstar.min        = zstarAnnual(:,5);
uvp.zstar.gof        = zstarMonthly_gof;

%% Save metrics

save(fullfile('.','data','processed','bcpmetrics_uvp.mat'),'uvp*')

end