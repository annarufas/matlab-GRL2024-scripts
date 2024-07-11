function [qSstMonthly,qChlaMonthly,qPar0Monthly] =... 
    extractHensonForcingDataForChosenLocations(isGlobalCalculation,...
    listLocalLats,listLocalLons)
   
% EXTRACTHENSONFORCINGDATAFORCHOSENLOCATIONS Extracts SST, chla and PAR0 
% for the locations defined by listLocalLats and listLocalLons.
%
%   INPUT: 
%       isGlobalCalculation - choice (1 or 0)
%       listLocalLats       - list of local latitudes
%       listLocalLons       - list of local longitudes
%
%   OUTPUT:
%       qSstMonthly  - monthly fields of SST (ºC)
%       qChlaMonthly - monthly fields of chla (mg m-3)
%       qPar0Monthly - monthly fields of PAR0 (W m-2)
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

    % Regrid SST, chla and PAR0 from the default grids used by AVHRR 
    % Pathfinder and SeaWiFS to the grid defined by qLats and qLons. 
    % Regridding is the process of interpolating from one grid resolution 
    % to a different one. 
    [Fchla,Fsst,Fpar0] = createInterpolants();
    
    [qX, qY, qT] = ndgrid(qLats, qLons, (1:12)'); 
    qSstMonthly = Fsst(qX, qY, qT);
    qSstMonthly = flipud(rot90(qSstMonthly)); % rotate the SST array so that it has the same lat/lon arrangement as chla and par0
    [qX, qY, qT] = ndgrid(qLons, qLats, (1:12)'); 
    qChlaMonthly = Fchla(qX, qY, qT);
    qPar0Monthly = Fpar0(qX, qY, qT);
    %figure(); pcolor(flipud(~rot90(qSstMonthly(:,:,1)))); caxis([-2 30]); shading interp;
   
    save(fullfile('.','data','processed','HensonForcingDataGlobal.mat'),...
        'qSstMonthly','qChlaMonthly','qPar0Monthly')

else
    
    qLats = listLocalLats;
    qLons = listLocalLons;
    nLocs = length(qLats);
    
    % Initiate output arrays. The size of each dimension is such that it
    % works with the Carr2002algorithm to calculate NPP
    qSstMonthly = NaN(nLocs,1,12);
    qChlaMonthly = NaN(nLocs,1,12);
    qPar0Monthly = NaN(nLocs,1,12);
    
    % Extract data for our chosen locations defined by qLats and qLons
    [Fchla,Fsst,Fpar0] = createInterpolants();
    for iLoc = 1:nLocs
        [qX, qY, qT] = ndgrid(qLats(iLoc), qLons(iLoc), (1:12)'); % query points for interpolation
        qSstMonthly(iLoc,:,:) = Fsst(qX, qY, qT);
        [qX, qY, qT] = ndgrid(qLons(iLoc), qLats(iLoc), (1:12)'); % REVERSE query points for interpolation
        qChlaMonthly(iLoc,:,:) = Fchla(qX, qY, qT);
        qPar0Monthly(iLoc,:,:) = Fpar0(qX, qY, qT);
    end

    save(fullfile('.','data','processed','HensonForcingDataLocal.mat'),...
        'qSstMonthly','qChlaMonthly','qPar0Monthly')

end

%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

% *************************************************************************
   
function [Fchla,Fsst,Fpar0] = createInterpolants()

    % SST from AVHRR Pathfinder v5.0 global monthly climatology, 1985–2001 
    load(fullfile('.','data','raw','sst_global.mat'),'sst','sst_lat','sst_lon'); % ºC

    % chla from SeaWiFS monthly climatology, 1997–2010 
    load(fullfile('.','data','raw','chla_seawifs_global.mat'),'chla','chla_lat','chla_lon') % mg m-3 (=ug L-1)

    % PAR0 from SeaWiFS monthly climatology, 1997–2010 
    load(fullfile('.','data','raw','par0_seawifs_global.mat'),'par0','par0_lat','par0_lon') % W m-2
    
    % Data grid
    [Xs, Ys, Ts] = ndgrid(sst_lat, sst_lon, (1:12)');
    [Xc, Yc, Tc] = ndgrid(chla_lon, chla_lat, (1:12)');
    [Xp, Yp, Tp] = ndgrid(par0_lon, par0_lat, (1:12)');
    
    % Interpolant -use first-order (linear) interpolation and extrapolation
    Fchla = griddedInterpolant(Xc, Yc, Tc, chla, 'linear'); 
    Fsst  = griddedInterpolant(Xs, Ys, Ts, sst, 'linear'); 
    Fpar0 = griddedInterpolant(Xp, Yp, Tp, par0, 'linear');

end % createInterpolants

end
