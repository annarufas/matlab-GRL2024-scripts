function [qTempMonthly] = extractMarsayForcingDataForChosenLocations(...
    isGlobalCalculation,listLocalLats,listLocalLons)
   
% EXTRACTMARSAYFORCINGDATAFORCHOSENLOCATIONS Extracts median temperature in
% the upper 500 m of the water column for the locations defined by 
% listLocalLats and listLocalLons.
%
%   INPUT: 
%       isGlobalCalculation - choice (1 or 0)
%       listLocalLats       - list of local latitudes
%       listLocalLons       - list of local longitudes
%
%   OUTPUT:
%       qTempMonthly  - monthly fields of temperature (ºC)
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

if isGlobalCalculation
    
    qLats = (-89.5:1:89.5)';
    qLons = (-179.5:1:179.5)';

    % Regrid temperature from the default grid used by WOA to the grid 
    % defined by qLats and qLons. Regridding is the process of interpolating 
    % from one grid resolution to a different one. 
    [Ftemp] = createInterpolants();
    [qX, qY, qT] = ndgrid(qLons, qLats, (1:12)); % query points for interpolation
    qTempMonthly = Ftemp(qX, qY, qT);
    
    save(fullfile('.','data','processed','MarsayForcingDataGlobal.mat'),'qTempMonthly')

else

    qLats = listLocalLats;
    qLons = listLocalLons;
    nLocs = length(qLats);
    
    qTempMonthly = NaN(nLocs,12);
    
    % Extract data for our chosen locations defined by qLats and qLons
    [Ftemp] = createInterpolants();
    for iLoc = 1:nLocs
        [qX, qY, qT] = ndgrid(qLons(iLoc), qLats(iLoc), (1:12)'); % query points for interpolation
        qTempMonthly(iLoc,:) = Ftemp(qX, qY, qT);
        disp(Ftemp(qX, qY, qT))
    end

    save(fullfile('.','data','processed','MarsayForcingDataLocal.mat'),'qTempMonthly')

end

%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

% *************************************************************************
   
function [Ftemp] = createInterpolants()
    
    % Climatological monthly mean temperature
    load(fullfile('.','data','raw','woa18_global.mat'),'WOA18temp','WOA18_lat','WOA18_lon','WOA18_zTemp') % ºC

    % Extract temperature in the first 500 m of the water column
    iDepth500 = find(WOA18_zTemp == 500);
    tempMonthlyUpp500 = squeeze(WOA18temp(:,:,(1:iDepth500),:));

    % Calculate the median temperature in the first 500 m of the water column
    tempMonthlyMedianUpp500 = squeeze(median(tempMonthlyUpp500(:,:,:,:),3,'omitnan'));

    % Data grid
    [X, Y, T] = ndgrid(WOA18_lon, WOA18_lat, (1:12));
    
    % Interpolant -use first-order (linear) interpolation and extrapolation
    Ftemp = griddedInterpolant(X, Y, T,tempMonthlyMedianUpp500, 'linear'); 

end % createInterpolants

end
