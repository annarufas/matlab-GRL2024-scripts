
% ======================================================================= %
%                                                                         %
% This script reads in particle concentration data from the UVP5          %
% instrument downloaded from the EcoTaxa repository and calculates POC    %
% flux using the method of Bisson et al. (2022). The script has 7         %
% sections:                                                               %
%   Section 1 - Presets.                                                  %
%   Section 2 - Read in EcoTaxa particle files and save the data into     %
%               .mat arrays.                                              %
%   Section 3 - Calculate POC flux.                                       %
%   Section 4 - Sort data into a unique array sorted by casts, depths,    %
%               months and location.                                      %
%   Section 5 - Calculate monthly averages and propagate error.           %
%   Section 6 - Calculate annual averages and propagate error.            % 
%   Section 7 – Save the data.                                            %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   WITH CODES PROVIDED BY K. BISSON, OREGON STATE                        %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 4 Jun 2024                                    %
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
isPocFluxReady = 0;
isAllDepths = 1;

if isAllDepths
    suffixFluxFilename = '_flux_data_all_depths.mat';
    ECOTAXA_VERTICAL_STEP = 5; % 5 m
    targetDepths = (2.5:ECOTAXA_VERTICAL_STEP:2005)'; % m
    filenameUvpProcessedDataset45sc = 'pocflux_bisson_45sc_monthly_and_annual_all_depths.mat';
else % use Bisson et al. (2022) depth definitions
    suffixFluxFilename = '_flux_data.mat';
    targetDepths = [7.5, 22.5, 47.5, 97.5, 147.5, 222.5, 297.5, 497.5, 997.5]; % m
    filenameUvpProcessedDataset45sc = 'pocflux_bisson_45sc_monthly_and_annual.mat';
end
nTargetDepths = numel(targetDepths);

% Particle size parameter declarations
NUM_SIZE_CLASSES = 45;
SIZE_STEP_PROGRESSION = 2^(1/3); % 2^(1/3) ~= 1.26 um

% Enter the coordinates that we have used to define our locations in the
% EcoTaxa's website map
NUM_LOCS = 6;
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
LAT_UPPER(4) = 32;     LAT_UPPER(5) = 23;      LAT_UPPER(6) = 80;
LAT_LOWER(4) = 29.5;   LAT_LOWER(5) = 22;      LAT_LOWER(6) = 78;
LON_RIGHT(4) = -62;    LON_RIGHT(5) = -157.5;  LON_RIGHT(6) = 5.5;
LON_LEFT(4) = -65;     LON_LEFT(5) = -158.5;   LON_LEFT(6) = 3.5;

% EcoTaxa folder definitions
SUFFIX_ECOTAXA_FOLDER_NAME = {'EqPac','OSP','PAPSO','BATSOFP','HOTALOHA','HAUSGARTEN'};
PREFIX_ECOTAXA_FOLDER_NAME = 'export_detailed_'; % common name part to all folders with UVP data

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - READ IN ECOTAXA PARTICLE FILES AND SAVE THE DATA INTO .MAT ARRAYS
% -------------------------------------------------------------------------

% Read in the *PART* files downloaded from the EcoTaxa website
if ~isEcotaxaDataReady
    readEcoTaxaParticleFiles(SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS)
end

% EcoTaxa size class definitions (um)
binEdges = zeros(NUM_SIZE_CLASSES+1,1);
binEdges(1) = 1; % 1 um
for i = 2:NUM_SIZE_CLASSES+1
    binEdges(i) = binEdges(i - 1) * SIZE_STEP_PROGRESSION; % 2^(1/3) ~= 1.26 um
end

binMiddle = zeros(NUM_SIZE_CLASSES,1);
for i = 1:NUM_SIZE_CLASSES
%     esdMiddle(i) = geomean(esdEdges(i:i+1));
    binMiddle(i) = (binEdges(i+1) + binEdges(i))./ 2;
end

binWidth = zeros(NUM_SIZE_CLASSES,1);
for i = 1:NUM_SIZE_CLASSES      
    binWidth(i) = binEdges(i+1) - binEdges(i);
end

% Transform units from um to mm
binEdges = binEdges*1e-3;
binMiddle = binMiddle*1e-3;
binWidth = binWidth*1e-3;

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CALCULATE POC FLUX
% -------------------------------------------------------------------------

% For all depths, this loop takes ~1 h

if ~isPocFluxReady

    % Number of times the modelled PSD flux will be sampled to calculate
    % uncertainty boundaries for POC flux
    nIterations = 100;

    tic

    for iLoc = 1:NUM_LOCS

        % The EcoTaxa data set
        load(fullfile('.','data','processed','UVP5',...
            strcat(SUFFIX_ECOTAXA_FOLDER_NAME{iLoc},'_particle_concentration.mat')))
        ET = pnumData;
        
        % Extract data only for the depths of interest
        matchingDepthRows = ismember(ET.Depth_m_, targetDepths);
        ET = ET(matchingDepthRows,:);

        % Crop the particle and sampled volume
        observedVol = ET.SampledVolume_L_; % L

        % Crop particle number
        observedNbl = zeros(height(ET),NUM_SIZE_CLASSES);
        for i = 1:NUM_SIZE_CLASSES
            observedNbl(:,i) = ET.(sizeClassLabels{i}); % # part. L-1
        end
        observedNbl(isnan(observedNbl)) = 0;

        for i = 1:length(observedNbl)
            
            % Only proceed if there are particles
            if (sum(observedNbl(i,:)) > 0)    
            
                % Calculate observed POC flux
                flux = calculatePocFlux(observedNbl(i,:)',binMiddle);
                ET.observedPocFlux(i) = flux; % mg C m-2 d-1
                
                % Calculate modelled POC flux
                observedN = round(observedNbl(i,:)'.*observedVol(i)); % # part. (integer)
                [estimNbl,C1,alpha,lambda] = calculateModelledPsd(observedN,...
                    observedVol(i),binWidth,binMiddle);
                flux = calculatePocFlux(estimNbl,binMiddle);
                
                ET.modelledPocFlux(i) = flux; % mg C m-2 d-1
                ET.C1(i)     = C1;
                ET.alpha(i)  = alpha;
                ET.lambda(i) = lambda;

                % Calculate POC flux for simulated PSD (this is to 
                % calculate POC flux uncertainty boundaries)
                simulatedNbl = generateRandomSamplesOfParticleNumber(...
                    observedNbl(i,:)',observedVol(i),nIterations,NUM_SIZE_CLASSES);

                fluxIter = NaN(nIterations,1);
                for iIter = 1:nIterations
                    if (nansum(simulatedNbl(iIter,:)) > 0)
                        fluxIter(iIter) = calculatePocFlux(simulatedNbl(iIter,:)',binMiddle);
                    end
                end

                fluxErr = std(fluxIter(:),'omitnan');
                ET.fluxErr(i) = fluxErr;

            end % if observedNbl > 0

        end

        % Save
        save(fullfile('.','data','processed','UVP5',...
            strcat(SUFFIX_ECOTAXA_FOLDER_NAME{iLoc},suffixFluxFilename)),'ET')

    end % iLoc

    toc
    
end % ~isPocFluxReady

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - SORT DATA INTO A UNIQUE ARRAY SORTED BY CASTS, DEPTHS, MONTHS
% AND LOCATION 
% -------------------------------------------------------------------------
    
% Calculate the number of casts (vertical sampling events) for each 
% location and month
castMonthlyDistrib = calculateNumberOfCasts(SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS,suffixFluxFilename);
maxNumCastsPerMonth = max(castMonthlyDistrib,[],'all');

% Combine local arrays of particle number concentrations into a single
% array with dimensions of cast x depth x month x location
[uvpFluxByCastAvg,uvpFluxByCastErr] = arrangePocFluxDataByLocMonthDepthCast(...
    SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS,maxNumCastsPerMonth,targetDepths,...
    nTargetDepths,suffixFluxFilename);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - CALCULATE MONTHLY AVERAGES AND PROPAGATE ERROR
% -------------------------------------------------------------------------

% Last dimension: 1=mean, 2=error, 3=nSamples
uvpMonthlyFlux = NaN(nTargetDepths,12,NUM_LOCS,3);

for iLoc = 1:NUM_LOCS
    
    for iMonth = 1:12
        nCasts = castMonthlyDistrib(iMonth,iLoc);
        
        if (nCasts > 0)
            for iDepth = 1:nTargetDepths

                % Number of casts that have reached that depth
                nCastsInDepth = sum(~isnan(uvpFluxByCastAvg(1:nCasts,iDepth,iMonth,iLoc)));
                
                if (nCastsInDepth > 0)
                    % Store the number of casts (or samples) at that depth
                    uvpMonthlyFlux(iDepth,iMonth,iLoc,3) = nCastsInDepth;
                
                    % The monthly POC flux value at that depth is calculated as 
                    % an average of the individual casts for that depth
                    uvpMonthlyFlux(iDepth,iMonth,iLoc,1) = mean(uvpFluxByCastAvg(1:nCasts,iDepth,iMonth,iLoc),'omitnan');

                    % Propagate errors from all the individual POC flux casts in 
                    % a month to monthly average of POC flux using "worstcase"
                    vals = uvpFluxByCastAvg(1:nCasts,iDepth,iMonth,iLoc);
                    err = uvpFluxByCastErr(1:nCasts,iDepth,iMonth,iLoc);

                    % Empty positions with NaN before propagating error
                    vals(isnan(vals)) = [];
                    err(isnan(err)) = [];
                    
                    if isequal(size(vals), size(err))
                        [~,~,f_LB,f_MID,f_UB,~,~] = worstcase(@(x) mean(x),vals,err);
                        uvpMonthlyFlux(iDepth,iMonth,iLoc,2) = f_UB - f_MID;
                    end
                    
                end % if nCastsInDepth > 0
            end % iDepth
        end % if nCasts > 0
    end % iMonth
end % iLoc

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 6 - CALCULATE ANNUAL AVERAGES AND PROPAGATE ERROR
% -------------------------------------------------------------------------

% Last dimension: 1=mean, 2=error, 3=max, 4=min
uvpAnnualFlux = NaN(nTargetDepths,NUM_LOCS,4); 

for iLoc = 1:NUM_LOCS
    for iDepth = 1:nTargetDepths
        nCastsInDepthForAllMonths = uvpMonthlyFlux(iDepth,:,iLoc,3);
            
        if (sum(isnan(nCastsInDepthForAllMonths)) > 0)
            vals = squeeze(uvpMonthlyFlux(iDepth,:,iLoc,1));
            err = squeeze(uvpMonthlyFlux(iDepth,:,iLoc,2));

            % Empty positions with NaN before propagating error
            vals(isnan(vals)) = [];
            err(isnan(err)) = [];
            nCastsInDepthForAllMonths(isnan(nCastsInDepthForAllMonths)) = [];
            
            if isequal(size(vals), size(err), size(nCastsInDepthForAllMonths))
                % Weighted mean (mw = ((mA*nA)+(mB*nB)+(mC*nC))/(nA+nB+nC))
                uvpAnnualFlux(iDepth,iLoc,1) = calculateWeightedAverage(vals,nCastsInDepthForAllMonths);

                % Error propagation
                [~,~,f_LB,f_MID,f_UB,~,~] = worstcase(@(x) calculateWeightedAverage(x,nCastsInDepthForAllMonths), vals', err');
                uvpAnnualFlux(iDepth,iLoc,2) = f_UB - f_MID;

                % Max and min values
                uvpAnnualFlux(iDepth,iLoc,3) = max(uvpMonthlyFlux(iDepth,:,iLoc,1));
                uvpAnnualFlux(iDepth,iLoc,4) = min(uvpMonthlyFlux(iDepth,:,iLoc,1));
            end

        end % if nCastsInDepthForAllMonths > 0
    end % iDepth
end % iLoc

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 7 - SAVE THE DATA
% -------------------------------------------------------------------------

% Save monthly and annual averages
save(fullfile('.','data','processed','UVP5',filenameUvpProcessedDataset45sc),...
    'uvpFluxByCastAvg','uvpFluxByCastErr','uvpMonthlyFlux','uvpAnnualFlux',...
    'targetDepths','castMonthlyDistrib','-v7.3')

fprintf('\nThe UVP5-derived POC flux dataset has been saved correctly.\n')

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

% *************************************************************************

function simulatedNbl = generateRandomSamplesOfParticleNumber(observedNbl,...
    observedVol,nIterations,NUM_SIZE_CLASSES)

% Assume observedNbl follows a negative binomial distribution and draw nIterations 
% samples from it to obtain simulatedN

simulatedN = NaN(nIterations,NUM_SIZE_CLASSES); 
simulatedNbl = NaN(nIterations,NUM_SIZE_CLASSES); 
for iSc = 1:NUM_SIZE_CLASSES
    if (observedNbl(iSc) > 0)
        simulatedN(:,iSc) = nbinrnd(observedNbl(iSc)+1/2,1/(observedVol+1),1,nIterations)'; % # part
        simulatedNbl(:,iSc) = simulatedN(:,iSc)./observedVol; % # part L-1
    end
end


% figure()
% plot(binMiddle(iFirst:iLast),simulatedPsd(1:10,:)')
% hold on
% plot(binMiddle(iFirst:iLast),observedNbl.*observedVol,'linewidth',3)
% set(gca,'yscale','log','xscale','log')
% xlabel('ESD (mm)')
% ylabel('Particle number')

end % generateRandomSamplesOfParticleNumber
                    
% *************************************************************************

function pocFlux = calculatePocFlux(psd,binMiddle) 

% Input:
%   psd, # part L-1
%   binWidth, mm
%   binMiddle, mm
% 
% Output:
%   pocFlux, mg C m-2 d-1

particleRadius = binMiddle./((4/3)*pi);
particleRadius = particleRadius.^(1/3);
particleEsd = particleRadius.*2; 

% Flux factor function from Eq. 5 in Bisson et al. (2022)
fluxFactor = @(d) 2.8.*d.^2.24;

pocFlux = nansum(psd.*fluxFactor(particleEsd));

end % calculatePocFlux

% *************************************************************************

function [estimNbl,C1,alpha,lambda] = calculateModelledPsd(observedN,observedVol,...
    binWidth,binMiddle) 

% Modelled particle number (# part./um), Eq. 1 in Bisson et al. (2022)
modelledN = @(x,d) x(1).*d(:).^-x(2).* exp(-d(:)./x(3)); % 'd' is esdMiddle and 'x' is an array containing C1, alpha and gamma

% log10 transform observedN and divide by bin width  
observedN_log10 = log10(observedN./binWidth); 

% Only consider populated size classes and remove the ones that are not
bad = find(isinf(observedN_log10)==1); 
observedN_log10(bad) = [];

binMiddle_good = binMiddle;
binMiddle_good(bad) = []; 

% Initial guess and bounds for model parameters
%       C1    alpha lambda
x0   = [40;   4;    0.1]; % # part L-1 / unitless / mm
xLow = [-inf; 0;    min(binMiddle_good)];
xUpp = [inf;  6;    max(binMiddle_good)]; 

% Optimise model parameters
xOptim = fminsearchbnd(@(x) norm(observedN_log10 - log10(modelledN(x,binMiddle_good))),... 
    x0, xLow, xUpp); 

C1     = xOptim(1);
alpha  = xOptim(2);
lambda = xOptim(3);

estimNbl = (modelledN(xOptim,binMiddle)./observedVol).*binWidth;

% Remove particle size classes that are not in the observed set
estimNbl(bad) = NaN;           

end % calculateModelledPsd

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

function listLocalEcoTaxaDirs = getLocalEcoTaxaDirectoryNames(...
    SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS)

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
    NUM_LOCS,suffixFluxFilename)

castMonthlyDistrib = zeros(12,NUM_LOCS);

for iLoc = 1:NUM_LOCS
 
    % The EcoTaxa data set
    load(fullfile('.','data','processed','UVP5',...
        strcat(SUFFIX_ECOTAXA_FOLDER_NAME{iLoc},suffixFluxFilename)),'ET')
    
    % Convert datetime columns and add 'month' and 'year' columns
    ET.date = datetime(ET.yyyy_mm_ddHh_mm,'format','yyyy-MM-dd');
    ET.dateString = datestr(ET.date, 'yyyy-mm-dd HH:MM:SS');

    for iMonth = 1:12
        monthFilter = month(ET.date) == iMonth;
        if (sum(monthFilter) > 0)
            % Get unique combinations of latitude, longitude and time
            [uniqueCombinations, ~, ~] = unique(ET{monthFilter,... 
                {'Latitude', 'Longitude','dateString'}}, 'rows');
            castMonthlyDistrib(iMonth,iLoc) = size(uniqueCombinations, 1);
        end
    end
    
end % iLoc

end % calculateNumberOfCasts

% *************************************************************************

function [uvpFluxByCastAvg,uvpFluxByCastErr] = arrangePocFluxDataByLocMonthDepthCast(...
    SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS,maxNumCastsPerMonth,targetDepths,...
    nTargetDepths,suffixFluxFilename)

uvpFluxByCastAvg = NaN(maxNumCastsPerMonth,nTargetDepths,12,NUM_LOCS);
uvpFluxByCastErr = NaN(maxNumCastsPerMonth,nTargetDepths,12,NUM_LOCS);
% uvpFluxByCastC1 = NaN(maxNumCastsPerMonth,nTargetDepths,12,NUM_LOCS);
% uvpFluxByCastAlpha = NaN(maxNumCastsPerMonth,nTargetDepths,12,NUM_LOCS);
% uvpFluxByCastGamma = NaN(maxNumCastsPerMonth,nTargetDepths,12,NUM_LOCS);

for iLoc = 1:NUM_LOCS
 
    % The EcoTaxa data set
    load(fullfile('.','data','processed','UVP5',...
        strcat(SUFFIX_ECOTAXA_FOLDER_NAME{iLoc},suffixFluxFilename)),'ET')
    
    % Convert datetime columns and add 'month' and 'year' columns
    ET.date = datetime(ET.yyyy_mm_ddHh_mm,'format','yyyy-MM-dd');
    ET.dateString = datestr(ET.date, 'yyyy-mm-dd HH:MM:SS');
    ET.month = month(ET.date);
    ET.year = year(ET.date);

    % Classify particle data by month, depth and cast
    for iMonth = 1:12
        monthFilter = ET.month == iMonth;

        if (sum(monthFilter) > 0)
            for iDepth = 1:nTargetDepths
                depthFilter = ET.Depth_m_ == targetDepths(iDepth) & monthFilter;

                if (sum(depthFilter) > 0)
                    
                    nInstances = numel(ET.modelledPocFlux(depthFilter));

                    % Store the filtered data in the output arrays
                    uvpFluxByCastAvg((1:nInstances),iDepth,iMonth,iLoc) = ET.modelledPocFlux(depthFilter);
                    uvpFluxByCastErr((1:nInstances),iDepth,iMonth,iLoc) = ET.fluxErr(depthFilter);

                end
                
            end % iDepth
        end
    end % iMonth
end % iLoc

end % arrangePocFluxDataByLocMonthDepthCast

% *************************************************************************
