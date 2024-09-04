
% ======================================================================= %
%                                                                         %
% This script reads in particle concentration data from the UVP5          %
% instrument downloaded from the EcoTaxa repository and calculates POC    %
% flux using the method of Bisson et al. (2022). The script has 3         %
% sections:                                                               %
%   Section 1 - Presets.                                                  %
%   Section 2 - Read in EcoTaxa particle files, extract relevant          % 
%               information for data processing (depths, ESD classes,     % 
%               casts) and sort data into a Matlab array.                 %
%   Section 3 - Calculate POC flux.                                       %           
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   WITH CODES PROVIDED BY K. BISSON, OREGON STATE                        %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 6 Jun 2024                                    %
%                                                                         %
% ======================================================================= %

close all; clear all; clc
addpath(genpath('./data/raw/'));
addpath(genpath('./data/processed/'));
addpath(genpath('./code/'));
addpath(genpath('./resources/external/'));
addpath(genpath('./resources/internal/'));

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Script options
isEcotaxaDataReady = 1;

% Parameter declarations
NUM_ESD_CLASSES = 45;
stepEsdProgression = 2^(1/3); % 2^(1/3) ~= 1.26 um
NUM_LOCS = 6;
ECOTAXA_VERTICAL_STEP = 5; % 5 m

% Enter the coordinates that we have used to define our locations in the
% EcoTaxa's website map
LAT_UPPER = zeros(NUM_LOCS,1);
LAT_UPPER = zeros(NUM_LOCS,1);
LON_RIGHT = zeros(NUM_LOCS,1);
LON_LEFT = zeros(NUM_LOCS,1);

% EqPac                % OSP                   % PAP-SO               
LAT_UPPER(1) = 4;      LAT_UPPER(2) = 51.5;    LAT_UPPER(3) = 49.5;   
LAT_LOWER(1) = -4;     LAT_LOWER(2) = 49.5;    LAT_LOWER(3) = 48.5;   
LON_RIGHT(1) = -148;   LON_RIGHT(2) = -144;    LON_RIGHT(3) = -16;    
LON_LEFT(1) = -152;    LON_LEFT(2) = -146;     LON_LEFT(3) = -17;     

% BATS/OFP             % HOT/ALOHA             % HAUSGARTEN  
LAT_UPPER(4) = 32;     LAT_UPPER(5) = 23 ;     LAT_UPPER(6) = 80;
LAT_LOWER(4) = 29.5;   LAT_LOWER(5) = 22;      LAT_LOWER(6) = 78;
LON_RIGHT(4) = -62;    LON_RIGHT(5) = -157.5;  LON_RIGHT(6) = 5.5;
LON_LEFT(4) = -65;     LON_LEFT(5) = -158.5;   LON_LEFT(6) = 3.5;

% EcoTaxa folder definitions
SUFFIX_ECOTAXA_FOLDER_NAME = {'EqPac','OSP','PAPSO','BATSOFP','HOTALOHA','HAUSGARTEN'};
PREFIX_ECOTAXA_FOLDER_NAME = 'export_detailed_'; % common nam part to all folders with UVP data

%%
% -------------------------------------------------------------------------
% SECTION 2 - READ IN ECOTAXA PARTICLE FILES, EXTRACT RELEVANT INFORMATION
% FOR DATA PROCESSING (DEPTHS, ESD CLASSES, CASTS) AND SORT DATA INTO A
% MATLAB ARRAY
% -------------------------------------------------------------------------

% Read in the *PART* files downloaded from the EcoTaxa website
if ~isEcotaxaDataReady
    readEcoTaxaParticleFiles(SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS)
end

% EcoTaxa depth definitions
ecotaxaDepths = (2.5:ECOTAXA_VERTICAL_STEP:2005)'; % m
nEcoTaxaDepths = numel(ecotaxaDepths);

% Bisson et al. (2022) target depths (p. 757)
bissonTargetDepths = [7.5, 22.5, 47.5, 97.5, 147.5, 222.5, 297.5, 497.5, 997.5]; % m
nBissonTargetDepths = numel(bissonTargetDepths);

% EcoTaxa 45 ESD class definitions
esdEdges = zeros(NUM_ESD_CLASSES+1,1);
esdEdges(1) = 1; % 1 um
for i = 2:NUM_ESD_CLASSES+1
    esdEdges(i) = esdEdges(i - 1) * stepEsdProgression; % 2^(1/3) ~= 1.26 um
end

esdMiddle = zeros(NUM_ESD_CLASSES,1);
for i = 1:NUM_ESD_CLASSES
    esdMiddle(i) = geomean(esdEdges(i:i+1));
end

esdBinWidth = zeros(NUM_ESD_CLASSES,1);
for i = 1:NUM_ESD_CLASSES      
    esdBinWidth(i) = esdEdges(i+1) - esdEdges(i);
end

% Calculate the number of casts (vertical sampling events) for each 
% location and month
castMonthlyDistrib = calculateNumberOfCasts(SUFFIX_ECOTAXA_FOLDER_NAME,...
                                            NUM_LOCS,...
                                            LAT_LOWER,LAT_UPPER,LON_LEFT,LON_RIGHT);
maxNumCasts = max(castMonthlyDistrib,[],'all');

% Combine local arrays of particle number concentrations into a single
% array with dimensions of ESD x cast x depth x month x location
[nbl,volsamp] = sortEcoTaxaDataByLocMonthDepthCast(SUFFIX_ECOTAXA_FOLDER_NAME,...
                                                   NUM_LOCS,...
                                                   NUM_ESD_CLASSES,...
                                                   maxNumCasts,...
                                                   ecotaxaDepths,...
                                                   nEcoTaxaDepths,...
                                                   LAT_LOWER,LAT_UPPER,LON_LEFT,LON_RIGHT);
                                               
% nbl: particle number, # part. L-1
% volsamp: sampled volume, L
                                               
%%
% -------------------------------------------------------------------------
% SECTION 3 - TRANSFORM UVP DATA INTO POC FLUX
% -------------------------------------------------------------------------

% Number of times the modelled POC flux will be sampled to calculate
% uncertainty boundaries
nIterations = 100;

% Initialise output arrays
C1             = NaN(maxNumCasts,nBissonTargetDepths,12,NUM_LOCS);
gamma          = NaN(maxNumCasts,nBissonTargetDepths,12,NUM_LOCS);
alpha          = NaN(maxNumCasts,nBissonTargetDepths,12,NUM_LOCS);
flux           = NaN(maxNumCasts,nBissonTargetDepths,12,NUM_LOCS); 
fluxIter       = NaN(nIterations,maxNumCasts,nBissonTargetDepths,12,NUM_LOCS);
fluxCastAvg    = NaN(maxNumCasts,nBissonTargetDepths,12,NUM_LOCS);
fluxCastErr    = NaN(maxNumCasts,nBissonTargetDepths,12,NUM_LOCS);

% The following loop takes around 2 h

tic

for iLoc = 1:NUM_LOCS
    
    for iMonth = 1:12
        nCasts = castMonthlyDistrib(iMonth,iLoc);
        
        if (nCasts > 0)
            for iDepth = 1:nBissonTargetDepths
                [~, idxBissonTargetDepth] = min(abs(ecotaxaDepths(:) - bissonTargetDepths(iDepth)));

                for iCast = 1:nCasts

                    observedNbl = nbl(:,iCast,idxBissonTargetDepth,iMonth,iLoc); % # part. L-1
                    observedNbl(isnan(observedNbl)) = 0;
                    
                    % Find the first and last positions occupied in the
                    % particle array
                    iFirst = find(observedNbl ~= 0, 1, 'first');
                    iLast = find(observedNbl ~= 0, 1, 'last');
                    nPopulatedEsdClasses = length(esdMiddle(iFirst:iLast));
                    
                    % Only proceed if there are particles in that cast at
                    % that depth
                    if (nPopulatedEsdClasses > 0)
                    
                        observedNbl = observedNbl(iFirst:iLast);                       % # part. L-1
                        observedVol = volsamp(iCast,idxBissonTargetDepth,iMonth,iLoc); % L
                        observedN   = round(observedNbl .* observedVol);               % # part. (integer)

                        % Calculate POC flux for observedN

                        if (sum(observedN(:)) > 0)        
                            [flux(iCast,iDepth,iMonth,iLoc),...
                            C1(iCast,iDepth,iMonth,iLoc),...
                            alpha(iCast,iDepth,iMonth,iLoc),...
                            gamma(iCast,iDepth,iMonth,iLoc)] =... 
                                calculatePocFlux(...
                                    observedN(:),esdBinWidth(iFirst:iLast),esdMiddle(iFirst:iLast));
                        end

                        % Calculate POC flux for simulatedN (this is to 
                        % calculate POC flux uncertainty boundaries)

                        simulatedN = generateRandomSamplesOfParticleNumber(...
                            observedNbl,observedVol,nIterations,nPopulatedEsdClasses);

                        for iIter = 1:nIterations
                            if (sum(simulatedN(iIter,:)) > 0)
                                [fluxIter(iIter,iCast,iDepth,iMonth,iLoc),~,~,~] =... 
                                    calculatePocFlux(...
                                        simulatedN(iIter,:)',esdBinWidth(iFirst:iLast),esdMiddle(iFirst:iLast));
                            end
                        end

                        fluxCastAvg(iCast,iDepth,iMonth,iLoc) = flux(iCast,iDepth,iMonth,iLoc);
                        fluxCastErr(iCast,iDepth,iMonth,iLoc) = std(fluxIter(:,iCast,iDepth,iMonth,iLoc),1,'omitnan');

                    end % nPopulatedEsdClasses > 0
                    
                end % iCast

%                 % The monthly POC flux value at that depth is calculated as 
%                 % an average of the individual casts for that depth
%                 uvpFluxMonthlyAvg(iDepth,iMonth,iLoc) = mean(fluxCastAvg(1:nCasts,iDepth,iMonth,iLoc),'omitnan');
%                 
%                 % Propagate errors from all the individual POC flux casts in 
%                 % a month to monthly average of POC flux using "worstcase"
% 
%                 thisDepthFluxCastsAvg = fluxCastAvg(1:nCasts,iDepth,iMonth,iLoc);
%                 thisDepthFluxCastsErr = fluxCastErr(1:nCasts,iDepth,iMonth,iLoc);
%                 
%                 thisDepthFluxCastsAvg(isnan(thisDepthFluxCastsAvg)) = 0;
%                 thisDepthFluxCastsErr(isnan(thisDepthFluxCastsErr)) = 0;
% 
%                 [x_LB,x_UB,f_LB,f_MID,f_UB,minus_percent,plus_percent] =... 
%                     worstcase(@(x) mean(x),thisDepthFluxCastsAvg,thisDepthFluxCastsErr);
%                 uvpFluxMonthlyErr(iDepth,iMonth,iLoc) = f_UB - f_MID;

            end % iDepth 
        end % if nCasts > 0
    end % iMonth
end % iLoc

toc

save(fullfile('.','data','processed','UVP5','Rufas_flux_45_size_classes.mat'),...
    'C1','gamma','alpha','flux','fluxIter','fluxCastAvg','fluxCastErr','-v7.3')


% Some checks
[maxC1,idxMaxC1] = max(C1,[],'all','linear');
[iC, jC, kC, lC] = ind2sub(size(C1), idxMaxC1);
[minAlpha,idxMinAlpha] = min(alpha,[],'all','linear');
[ia, ja, ka, la] = ind2sub(size(alpha), idxMinAlpha);

%%
% -------------------------------------------------------------------------
% SECTION 4 - CALCULATE MONTHLY AVERAGES AND PROPAGATE ERROR
% -------------------------------------------------------------------------

uvpMonthlyFluxAvg     = NaN(nBissonTargetDepths,12,NUM_LOCS);
uvpMonthlyFluxErr     = NaN(nBissonTargetDepths,12,NUM_LOCS);
uvpMonthlyFluxSamples = NaN(nBissonTargetDepths,12,NUM_LOCS);

for iLoc = 1:NUM_LOCS
    
    for iMonth = 1:12
        nCasts = castMonthlyDistrib(iMonth,iLoc);
        
        if (nCasts > 0)
            for iDepth = 1:nBissonTargetDepths

                % Number of casts that have reacehd that depth
                
                uvpMonthlyFluxSamples(iDepth,iMonth,iLoc) = sum(~isnan(...
                    fluxCastAvg(1:nCasts,iDepth,iMonth,iLoc)));
                
                % The monthly POC flux value at that depth is calculated as 
                % an average of the individual casts for that depth
                
                uvpMonthlyFluxAvg(iDepth,iMonth,iLoc) = mean(fluxCastAvg(1:nCasts,iDepth,iMonth,iLoc),'omitnan');
                
                % Propagate errors from all the individual POC flux casts in 
                % a month to monthly average of POC flux using "worstcase"

                vals = fluxCastAvg(1:nCasts,iDepth,iMonth,iLoc);
                err = fluxCastErr(1:nCasts,iDepth,iMonth,iLoc);
                
                % Empty positions with no values
                vals(isnan(vals)) = [];
                err(isnan(err)) = [];

                [x_LB,x_UB,f_LB,f_MID,f_UB,minus_percent,plus_percent] =... 
                    worstcase(@(x) mean(x),vals,err);
                uvpMonthlyFluxErr(iDepth,iMonth,iLoc) = f_UB - f_MID;

            end % iDepth
        end % if nCasts > 0
    end % iMonth
end % iLoc

%%
% -------------------------------------------------------------------------
% SECTION 5 - CALCULATE ANNUAL AVERAGES AND PROPAGATE ERROR
% -------------------------------------------------------------------------

% Last dimension: 1=mean, 2=total error, 3=max, 4=min
uvpAnnualFlux = NaN(nBissonTargetDepths,NUM_LOCS,4); 

for iLoc = 1:NUM_LOCS
    
    for iDepth = 1:nBissonTargetDepths
        nCastsInDepth = uvpMonthlyFluxSamples(iDepth,iMonth,iLoc);

        if (nCastsInDepth > 0)

            vals = squeeze(uvpMonthlyFluxAvg(iDepth,:,iLoc));
            err = squeeze(uvpMonthlyFluxErr(iDepth,:,iLoc));

            % Empty positions with no values
            vals(isnan(vals)) = [];
            err(isnan(err)) = [];

            % Weighted mean (mw = ((mA*nA)+(mB*nB)+(mC*nC))/(nA+nB+nC))
            uvpAnnualFlux(iDepth,iLoc,1) = calculateWeightedAverage(vals,nCastsInDepth);

            % Error propagation
            [x_LB, x_UB, f_LB, f_MID, f_UB, minus_percent, plus_percent] = ...
                worstcase(@(x) calculateWeightedAverage(x, nCastsInDepth), vals', err');
            uvpAnnualFlux(iDepth,iLoc,2) = f_UB - f_MID;

            % Max and min values
            uvpAnnualFlux(iDepth,iLoc,3) = max(uvpMonthlyFluxErr(iDepth,:,iLoc,1));
            uvpAnnualFlux(iDepth,iLoc,4) = min(uvpMonthlyFluxErr(iDepth,:,iLoc,1));

        end
    end % iDepth
end % iLoc

%%
% -------------------------------------------------------------------------
% SECTION 6 - FIND MATCHUPS BETWEEN THE UVP5-DERIVED ESTIMATES AND THE TRAP
% AND RADIONUCLIDE OBSERVATIONS
% -------------------------------------------------------------------------

% Extract observations from those months where there are UVP5 observations.
% Use the list of observed depths to match and extract UVP5 flux estimates.

matchObsFluxesMean   = NaN(MAX_NUM_DEPTHS_PER_PROFILE_COMPILATION,12,NUM_LOCS); 
matchObsFluxesErrTot = NaN(MAX_NUM_DEPTHS_PER_PROFILE_COMPILATION,12,NUM_LOCS); 
matchObsFluxesDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE_COMPILATION,12,NUM_LOCS);               

matchEstFluxesMean   = NaN(MAX_NUM_DEPTHS_PER_PROFILE_COMPILATION,12,NUM_LOCS); 
matchEstFluxesErrTot = NaN(MAX_NUM_DEPTHS_PER_PROFILE_COMPILATION,12,NUM_LOCS);
matchEstFluxesDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE_COMPILATION,12,NUM_LOCS);

for iLoc = 1:NUM_LOCS
    for iMonth = 1:12
        
        if (nCastsMonthly(iMonth,iLoc) > 0)
                         
            matchObsFluxesMean(:,iMonth,iLoc) = obsMonthlyProfileAvg(:,iMonth,iLoc);    
            matchObsFluxesErrTot(:,iMonth,iLoc) = obsMonthlyProfileErrTot(:,iMonth,iLoc);  
            matchObsFluxesDepths(:,iMonth,iLoc) = obsMonthlyProfileDepths(:,iMonth,iLoc);

            for iUniqueDepth = 1:nUniqueObsDepths(iMonth,iLoc)
                
                zo = matchObsFluxesDepths(iUniqueDepth,iMonth,iLoc);
                
                % Only proceed if the observed depth is smaller than the 
                % lower limit depth recorded by the UVP + 10 m
                if (zo <= (uvpDepthsPerCast(end) + 10)) 
                    [diff,idx] = min(abs(uvpDepthsPerCast(:)-zo));
                    if (diff <= MATCHUP_DEPTH_SEARCH)
                        matchEstFluxesDepths(iUniqueDepth,iMonth,iLoc) =...
                            uvpDepthsPerCast(idx);
                        matchEstFluxesMean(iUniqueDepth,iMonth,iLoc) =...
                            uvpMonthlyFlux(idx,iMonth,iLoc,1);
                        matchEstFluxesErrTot(iUniqueDepth,iMonth,iLoc) =... 
                            uvpMonthlyFlux(idx,iMonth,iLoc,4);
                    end
                end
                
            end % iUniqueDepth
        end % nCastsMonthly(iMonth,iLoc) > 0
    end % iMonth
end % iLoc

nMatchups = nnz(~isnan(matchEstFluxesMean));
fprintf('%\nd matchups were found between obserevd and estimated values of POC flux.\n', nMatchups)

%%
% -------------------------------------------------------------------------
% SECTION 7 - SAVE THE DATA
% -------------------------------------------------------------------------

% Save monthly and annual averages
save(fullfile('.','data','processed',filenameUvpProcessedDataset45sc),...
    'uvpFluxByCast','uvpMonthlyFlux','uvpAnnualFlux','uvpDepthsPerCast',...
    'nCastsMonthly','nDepthsUvp')

% Save matchup data
save(fullfile('.','data','processed',filenameObsVsEstMatchups),...
    'matchEstFluxesMean','matchEstFluxesErrTot','matchEstFluxesDepths',...
    'matchObsFluxesMean','matchObsFluxesErrTot','matchObsFluxesDepths',...
    'nMatchups')

fprintf('\nThe UVP5-derived POC flux dataset has been saved correctly.\n')


% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

% *************************************************************************

function simulatedN = generateRandomSamplesOfParticleNumber(observedNbl,...
    observedVol,nIterations,nPopulatedEsdClasses)

% Assume observedNbl follows a negative binomial distribution and draw nIterations 
% samples from it to obtain simulatedN

simulatedN  = NaN(nIterations,nPopulatedEsdClasses); 
for iEsd = 1:nPopulatedEsdClasses
    simulatedN(:,iEsd) = nbinrnd(observedNbl(iEsd)+1/2,1/(observedVol+1),1,nIterations)';
end

% figure()
% plot(esdMiddle(iFirst:iLast),simulatedN(1:10,:)')
% hold on
% plot(esdMiddle(iFirst:iLast),observedNbl.*observedVol,'linewidth',3)
% set(gca,'yscale','log','xscale','log')
% xlabel('ESD (um)')
% ylabel('Particle number')

end % generateRandomSamplesOfParticleNumber
                    
% *************************************************************************

function [pocFlux,C1,alpha,gamma] = calculatePocFlux(particleNumber,...
    binWidth,binMiddle) 

% Modelled particle number (# part./um), Eq. 1 in Bisson et al. (2022)
modelledN = @(x,d) x(1).*d(:).^-x(2).* exp(-d(:)./x(3)); % 'd' is esdMiddle in um and 'x' is an array containing C1, alpha and gamma

% Initial guess and bounds for model parameters
%       C1    alpha gamma
x0   = [40;   4;    100]; % # part L-1 / unitless / um
xLow = [-inf; 0;    min(binMiddle)];
xUpp = [inf;  6;    max(binMiddle)];                  

% Flux factor function from Eq. 5 in Bisson et al. (2022)
fluxFactor = @(d) 2.8.*(d(:).*1e-3).^2.24; % 'd' is esdMiddle in mm

% log10 transform observedN and divide by bin width  
observedN_log10 = log10(particleNumber./binWidth); 

% Only consider populated size classes and remove the ones that are not
bad = find(isinf(observedN_log10)==1); 
observedN_log10(bad) = [];

binMiddle_good = binMiddle;
binMiddle_good(bad) = []; 

xOptim = fminsearchbnd(@(x) norm(observedN_log10 - log10(modelledN(x,binMiddle_good))),... 
    x0, xLow, xUpp); % optimised parameters C1, alpha and gamma

C1    = xOptim(1);
alpha = xOptim(2);
gamma = xOptim(3);

estimPsd = modelledN(xOptim,binMiddle).*binWidth;

pocFlux = nansum(estimPsd.*fluxFactor(binMiddle)); % mg C m-2 d-1

end % calculatePocFlux
                             
% *************************************************************************                             

function readEcoTaxaParticleFiles(SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS)
    
for iLoc = 1:NUM_LOCS
 
    listLocalEcoTaxaDirs = getLocalEcoTaxaDirectoryNames(SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS);
    pathDir = fullfile('.','data','raw','UVP5',listLocalEcoTaxaDirs{iLoc});

    % Metadata
    metadataFile = dir(fullfile(pathDir,'*metadata*.tsv'));
    M = readtable(fullfile(pathDir,metadataFile.name),...
        'FileType','text','Delimiter','tab','TreatAsEmpty',{'N/A','n/a'}); % metadata table
    
    % Keep only relevant metadata columns
    M = M(:, {'profile', 'Cruise', 'Longitude', 'Latitude'});

    listOfParticleFiles = dir(fullfile(pathDir,'*PAR*.tsv')); % 'PAR' stands for 'particle'
    nSitesSampledForParticles = length(listOfParticleFiles);

    for iSite = 1:nSitesSampledForParticles

         thisParticleFile = listOfParticleFiles(iSite).name;
         T = readtable(fullfile(pathDir,thisParticleFile),...
             'FileType','text','Delimiter','tab','TreatAsEmpty',{'N/A','n/a'}); % file table

         % On the first iteration, identify particle size classes
         if (iSite == 1)
            idxsSizeClasses = strncmp(T.Properties.VariableNames,'LPM_',length('LPM_')); % 'LPM' stands for particle size class
            sizeClassLabels = T.Properties.VariableNames(idxsSizeClasses)'; 
         end

         % Crop relevant particle information
         Tt3 = T(:,idxsSizeClasses); % crop file table
         Tt2 = T(:,(3:5)); % date, depth and sampled volume
         Tt1 = repmat(M(iSite,:),[height(Tt3) 1]); % repeat profile, cruise, longitude and latitude data

         % Combine metadata, particle info, and particle size data
         P = [Tt1,Tt2,Tt3]; 

         if (iSite == 1 && ~isempty(P))
            vars = P.Properties.VariableNames;
            pnumData = cell2table(cell(0,length(vars)), 'VariableNames', vars);
            pnumData = P;
         elseif (iSite > 1 && ~isempty(P))
            % Append data from the current file to the existing particle data
            pnumData = [pnumData; P];
         end

    end % iLocSite
    
    save(fullfile('.','data','processed','UVP5',...
        strcat(SUFFIX_ECOTAXA_FOLDER_NAME{iLoc},'_particle_concentration.mat')),...
        'pnumData','sizeClassLabels','vars')

end % iLoc

end % readEcoTaxaParticleFiles

% *************************************************************************

function listLocalEcoTaxaDirs = getLocalEcoTaxaDirectoryNames(SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS)

% EcoTaxa directories are named using a specific structure consisting of a
% "prefix", "date of data download" and "suffix". This section manages
% data downloaded on various dates and organises the list of EcoTaxa
% directories according to the order specified by the parameter
% SUFFIX_ECOTAXA_FOLDER_NAME.

% Get a list of all directories in the search path and filter out non-directory entries
pathEcoTaxaDirs = dir(fullfile('.','data','raw','UVP5'));
pathEcoTaxaDirs = pathEcoTaxaDirs([pathEcoTaxaDirs.isdir]);
nameEcoTaxaDirs = {pathEcoTaxaDirs(3:end).name}; % remove '.' and '..' from the list of directories

% Create a mapping from each suffix to its index in the suffix array
suffixToIndexMap = containers.Map(SUFFIX_ECOTAXA_FOLDER_NAME, 1:NUM_LOCS);

% Initialise an array to store the reordered directory names
listLocalEcoTaxaDirs = cell(NUM_LOCS,1);
for i = 1:numel(nameEcoTaxaDirs)
    splitName = strsplit(nameEcoTaxaDirs{i}, '_');
    thisSubdirSuffix = splitName{4}; % get fourth part
    
    % Find the index of the suffix in the suffix array using the mapping
    desiredIndex = suffixToIndexMap(thisSubdirSuffix);
    
    % Place the subdirectory name in the desired position in the reordered array
    listLocalEcoTaxaDirs{desiredIndex} = nameEcoTaxaDirs{i};
end

end % getLocalEcoTaxaDirectoryNames

% *************************************************************************

function castMonthlyDistrib = calculateNumberOfCasts(SUFFIX_ECOTAXA_FOLDER_NAME,...
    NUM_LOCS,LAT_LOWER,LAT_UPPER,LON_LEFT,LON_RIGHT)

castMonthlyDistrib = zeros(12,NUM_LOCS);

for iLoc = 1:NUM_LOCS
 
    % The EcoTaxa data set
    load(fullfile('.','data','processed','UVP5',...
        strcat(SUFFIX_ECOTAXA_FOLDER_NAME{iLoc},'_particle_concentration.mat')))
    ET = pnumData;
    
    % Convert datetime columns and add 'month' and 'year' columns
    ET.date = datetime(ET.yyyy_mm_ddHh_mm,'format','yyyy-MM-dd');
    ET.dateString = datestr(ET.date, 'yyyy-mm-dd HH:MM:SS');

    % Define latitude and longitude bounds for the location
    latlonFilter = ET.Latitude >= LAT_LOWER(iLoc) & ET.Latitude < LAT_UPPER(iLoc) & ...
                   ET.Longitude >= LON_LEFT(iLoc) & ET.Longitude < LON_RIGHT(iLoc);

    if (sum(latlonFilter) > 0)
        for iMonth = 1:12
            monthFilter = month(ET.date) == iMonth & latlonFilter;
            if (sum(monthFilter) > 0)
                % Get unique combinations of latitude, longitude and time
                [uniqueCombinations, ~, ~] = unique(ET{monthFilter,... 
                    {'Latitude', 'Longitude','dateString'}}, 'rows');
                castMonthlyDistrib(iMonth,iLoc) = size(uniqueCombinations, 1);
            end
        end
    end
    
end % iLoc

end % calculateNumberOfCasts

% *************************************************************************

function [nbl,volsamp] = sortEcoTaxaDataByLocMonthDepthCast(...
    SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS,NUM_ESD_CLASSES,maxNumCasts,...
    ecotaxaDepths,nEcoTaxaDepths,LAT_LOWER,LAT_UPPER,LON_LEFT,LON_RIGHT)

nbl = NaN(NUM_ESD_CLASSES,maxNumCasts,nEcoTaxaDepths,12,NUM_LOCS); % # part. L-1
volsamp = NaN(maxNumCasts,nEcoTaxaDepths,12,NUM_LOCS); % L

for iLoc = 1:NUM_LOCS
 
    % The EcoTaxa data set
    load(fullfile('.','data','processed','UVP5',...
        strcat(SUFFIX_ECOTAXA_FOLDER_NAME{iLoc},'_particle_concentration.mat')))
    ET = pnumData;
    
    % Convert datetime columns and add 'month' and 'year' columns
    ET.date = datetime(ET.yyyy_mm_ddHh_mm,'format','yyyy-MM-dd');
    ET.dateString = datestr(ET.date, 'yyyy-mm-dd HH:MM:SS');
    ET.month = month(ET.date);
    ET.year = year(ET.date);

    % Crop the particle and sampled volume data
    ETpnum = zeros(height(ET),NUM_ESD_CLASSES);
    for i = 1:NUM_ESD_CLASSES
        ETpnum(:,i) = ET.(sizeClassLabels{i}); % # part. L-1
    end
    ETvol = ET.SampledVolume_L_; % L

    % Define latitude and longitude bounds for the location
    latlonFilter = ET.Latitude >= LAT_LOWER(iLoc) & ET.Latitude < LAT_UPPER(iLoc) & ...
                   ET.Longitude >= LON_LEFT(iLoc) & ET.Longitude < LON_RIGHT(iLoc);

    % Classify particle data by month, depth and cast
    if (sum(latlonFilter) > 0)
        
        for iMonth = 1:12
            monthFilter = ET.month == iMonth & latlonFilter;
            
            if (sum(monthFilter) > 0)

                for iDepth = 1:nEcoTaxaDepths
                    depthFilter = ET.Depth_m_ == ecotaxaDepths(iDepth) &... 
                                  monthFilter;

                    if (sum(depthFilter) > 0)
                        
                        % Get volume 
                        nInstances = numel(ETvol(depthFilter));
                        volsamp((1:nInstances),iDepth,iMonth,iLoc) = ETvol(depthFilter);
                        
                        for iEsd = 1:NUM_ESD_CLASSES
                            esdFilter = ~isnan(ET.(sizeClassLabels{iEsd})) &... 
                                        ET.(sizeClassLabels{iEsd}) > 0 &... 
                                        depthFilter;
     
                            if (sum(esdFilter) > 0)
                                
                                % Get particle number
                                nInstances = numel(ETpnum(esdFilter,iEsd));
                                nbl(iEsd,(1:nInstances),iDepth,iMonth,iLoc) = ETpnum(esdFilter,iEsd);
                                
                            end
                
                        end % iEsd
                    end
                end % iDepth
            end
        end % iMonth
    end
end % iLoc

end % sortEcoTaxaDataByLocMonthDepthCast

% *************************************************************************
