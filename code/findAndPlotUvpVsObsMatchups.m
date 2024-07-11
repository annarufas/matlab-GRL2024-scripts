% ======================================================================= %
%                                                                         %
% This script finds matchups between the POC flux data from the trap and  %
% radionuclide compilation and the POC flux estimated from the UVP5       % 
% dataset. The script is structured into 3 sections:                      %
%   Section 1 - Presets.                                                  %           
%   Section 2 - Find matchups between the UVP5-derived estimates and the  %
%               trap and radionuclide compialtion.                        %
%   Section 3 – Plot matchups figure (Figure 3).                          %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
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

filenameUvpProcessedDataset45sc = 'pocflux_bisson_45sc_monthly_and_annual_all_depths.mat';
filenameMonthlyPocFlux          = 'pocflux_compilation_monthly.mat';
filenameTimeseriesInformation   = 'timeseries_station_information.mat';

MAX_NUM_DEPTHS_PER_PROFILE_COMPILATION = 100; % in the compiled dataset
MATCHUP_DEPTH_SEARCH = 5; % m

% Load the UVP5 data
load(fullfile('.','data','processed','UVP5',filenameUvpProcessedDataset45sc),...
    'targetDepths','castMonthlyDistrib','uvpMonthlyFlux')

% Load the trap and radionuclide compilation
load(fullfile('.','data','processed',filenameMonthlyPocFlux),...
    'obsMonthlyProfileAvg','obsMonthlyProfileErrTot','obsMonthlyProfileDepths',...
    'nUniqueObsDepths')

% Load station information
load(fullfile('.','data','processed',filenameTimeseriesInformation))

%%
% -------------------------------------------------------------------------
% SECTION 2 - FIND MATCHUPS BETWEEN THE UVP5-DERIVED ESTIMATES AND THE TRAP
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
        
        if (castMonthlyDistrib(iMonth,iLoc) > 0)
                         
            matchObsFluxesMean(:,iMonth,iLoc) = obsMonthlyProfileAvg(:,iMonth,iLoc);    
            matchObsFluxesErrTot(:,iMonth,iLoc) = obsMonthlyProfileErrTot(:,iMonth,iLoc);  
            matchObsFluxesDepths(:,iMonth,iLoc) = obsMonthlyProfileDepths(:,iMonth,iLoc);

            for iUniqueDepth = 1:nUniqueObsDepths(iMonth,iLoc)
                
                zo = matchObsFluxesDepths(iUniqueDepth,iMonth,iLoc);
                
                % Only proceed if the observed depth is smaller than the 
                % lower limit depth recorded by the UVP + 10 m
                if (zo <= (targetDepths(end) + 10)) 
                    [diff,idx] = min(abs(targetDepths(:)-zo));
                    
                    if (diff <= MATCHUP_DEPTH_SEARCH)
                        matchEstFluxesDepths(iUniqueDepth,iMonth,iLoc) =...
                            targetDepths(idx);
                        matchEstFluxesMean(iUniqueDepth,iMonth,iLoc) =...
                            uvpMonthlyFlux(idx,iMonth,iLoc,1);
                        matchEstFluxesErrTot(iUniqueDepth,iMonth,iLoc) =... 
                            uvpMonthlyFlux(idx,iMonth,iLoc,2);
                    end
                end
                
            end % iUniqueDepth
        end % castMonthlyDistrib(iMonth,iLoc) > 0
    end % iMonth
end % iLoc

nMatchups = nnz(~isnan(matchEstFluxesMean));
fprintf('%d matchups were found between obserevd and estimated values of POC flux.\n', nMatchups)


obsmeas = reshape(obsMonthlyProfileAvg, [], 1);
monthlyNumebrOfEntriesObserved = sum(~isnan(obsmeas), 'all');

uvp = squeeze(uvpMonthlyFlux(:,:,:,1));
uvpresh = reshape(uvp, [], 1);
monthlyNumberOfEntriesEstimated = sum(~isnan(uvpresh), 'all');

%%
% -------------------------------------------------------------------------
% SECTION 3 - PLOT FIGURE 3 (MATCHUPS)
% -------------------------------------------------------------------------

% Compute overall Spearman rank correlation coefficient 
% Get rid of NaN so that the computation can proceed.
idxNonNaN = ~isnan(matchEstFluxesMean(:,:,:));
x = matchEstFluxesMean(:,:,:);
y = matchObsFluxesMean(:,:,:);
[rho,pval] = corr(x(idxNonNaN),y(idxNonNaN),'Type','Spearman');

% Compute confidence intervals using an OLS regression determined via 
% non-parametric bootstrapping (as in Fender et al. 2019)
% X = [];
% Y = [];
% for iStation = 1:NUM_LOCS
%     switch iStation
%         case 1
%             iLoc = 5; % HOT/ALOHA
%         case 2
%             iLoc = 4; % BATS/OFP
%         case 3
%             iLoc = 1; % EqPac
%         case 4
%             iLoc = 3; % PAP-SO
%         case 5
%             iLoc = 2; % OSP
%         case 6
%             iLoc = 6; % HAUSGARTEN    
%     end
%     x = squeeze(matchObsFluxesMean(:,:,iLoc));
%     y = squeeze(matchEstFluxesMean(:,:,iLoc));
%     idxNonNaN = ~isnan(matchEstFluxesMean(:,:,iLoc));
%     xm = x(idxNonNaN);
%     ym = y(idxNonNaN);
%     X = [X; xm];
%     Y = [Y; ym];
% end
% 
% % Perform non-parametric bootstrapping
% % yPredictedBootstrap = zeros(nBootstrapSamples, length(X));
% nBootstrapSamples = 1000;
% xValues = linspace(0, 150, 100); % same as in @mypredict
% yPredictedBootstrap = bootstrp(nBootstrapSamples, @mypredict, [X, Y]);
% lowerCI = prctile(yPredictedBootstrap, 2.5, 1);
% upperCI = prctile(yPredictedBootstrap, 97.5, 1);

% .........................................................................

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.28 0.36],'Color','w')
coloursLocs = brewermap(NUM_LOCS,'*Set1');
axh = axes('Position', [0.11 0.11 0.62 0.72]);
set(axh,'ColorOrder',coloursLocs,'NextPlot','replacechildren');

% Plot error bars first
for iStation = 1:NUM_LOCS

    % Re-order
    switch iStation
        case 1
            iLoc = 5; % HOT/ALOHA
        case 2
            iLoc = 4; % BATS/OFP
        case 3
            iLoc = 1; % EqPac
        case 4
            iLoc = 3; % PAP-SO
        case 5
            iLoc = 2; % OSP
        case 6
            iLoc = 6; % HAUSGARTEN  
    end 

    x = squeeze(matchObsFluxesMean(:,:,iLoc));
    y = squeeze(matchEstFluxesMean(:,:,iLoc));
    xerr = squeeze(matchObsFluxesErrTot(:,:,iLoc));
    yerr = squeeze(matchEstFluxesErrTot(:,:,iLoc));
    idxNonNaN = ~isnan(matchEstFluxesMean(:,:,iLoc));
    scatter(axh,NaN,NaN,1,'o','MarkerEdgeColor','k','MarkerFaceColor',coloursLocs(iLoc,:),...
        'LineWidth',0.5,'HandleVisibility','off');
    hold on
    
    eb(1) = errorbar(x(idxNonNaN),y(idxNonNaN),xerr(idxNonNaN),'horizontal',...
        'LineStyle','none','HandleVisibility','off');
    eb(2) = errorbar(x(idxNonNaN),y(idxNonNaN),yerr(idxNonNaN),'vertical',...
        'LineStyle','none','HandleVisibility','off');
    eb(1).CapSize = 0;
    eb(2).CapSize = 0;
    set(eb,'Color',coloursLocs(iLoc,:),'LineWidth',1)
    
    % Set transparency level (0:1)
    alpha = 0.3;  
    set([eb(1).Bar,eb(1).Line],'ColorType','truecoloralpha',...
        'ColorData',[eb(1).Line.ColorData(1:3); 255*alpha])
    set([eb(2).Bar,eb(2).Line],'ColorType','truecoloralpha',...
        'ColorData',[eb(2).Line.ColorData(1:3); 255*alpha])
    hold on

end
hold on

% Plot scatters on top
for iStation = 1:NUM_LOCS

    % Re-order
    switch iStation
        case 1
            iLoc = 5; % HOT/ALOHA
        case 2
            iLoc = 4; % BATS/OFP
        case 3
            iLoc = 1; % EqPac
        case 4
            iLoc = 3; % PAP-SO
        case 5
            iLoc = 2; % OSP
        case 6
            iLoc = 6; % HAUSGARTEN
    end 

    x = squeeze(matchObsFluxesMean(:,:,iLoc));
    y = squeeze(matchEstFluxesMean(:,:,iLoc));
    idxNonNaN = ~isnan(matchEstFluxesMean(:,:,iLoc));
    scatter(axh,x(idxNonNaN),y(idxNonNaN),50,'o','MarkerEdgeColor','k',...
        'MarkerFaceColor',coloursLocs(iLoc,:),'LineWidth',0.5);
    hold on

    % Plot 1:1 reference line
    if (iStation == NUM_LOCS) 
        hline = refline(axh,1,0); 
        hline.Color = 'k';
        hline.LineStyle = '--';
    end
    hold on

end
hold on
box on

% Plot 95% confidence intervals
% fill([xValues,fliplr(xValues)], [lowerCI,fliplr(upperCI)],...
%     'r','FaceAlpha',0.3,'EdgeColor','none');

% Calculate xlim(2) and ylim(2)
% maxObs = max(matchObsFluxesMean,[],'all','omitnan');
% maxEst = max(matchEstFluxesMean,[],'all','omitnan');
% maxFluxVal = max(maxObs,maxEst);   

ylim([0 230])
xlim([0 230])
% yticks([0 25 50 75 100 125 150]) 
% xticks([0 25 50 75 100 125 150]) 

% Add statistics
xt = max(xlim)-0.02*max(xlim); 
yt = max(ylim)-0.02*max(ylim);
text(yt,xt,... % text position relative to axis
    strcat('{\it r} =',{' '},num2str(rho,'%.2f'),{' '},'({\it p} =',{' '},num2str(pval,'%.2f'),')'),...
    'FontSize',12,'Horiz','right','Vert','top');
text(yt,xt-14,... % text position relative to axis
    strcat('{\it N} =',{' '},num2str(nMatchups,'%.0f')),...
    'FontSize',12,'Horiz','right','Vert','top');

% Add legend
lg = legend(axh,{'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN','1:1 line'});  
lg.Position(1) = 0.74; lg.Position(2) = 0.55;
lg.ItemTokenSize = [15,1];
lg.FontSize = 12;
set(lg,'Box','off')

ylabel('Estimated (UVP5)')
xlabel('Observed (sediment traps & radionuclides)')
title('POC flux (mg m^{-2} d^{-1})')

xh = get(axh,'xlabel'); 
px = get(xh,'position');
px(2) = px(2) - 2;       
set(xh,'position',px)

yh = get(axh,'ylabel'); 
py = get(yh,'position');
py(1) = py(1) - 3;
set(yh,'position',py)

th = get(axh,'title'); 
tp = get(th,'position');
tp(2) = tp(2) + 2;        
set(th,'position',tp)

set(axh,'FontSize',12)

exportgraphics(gcf,fullfile('.','figures','matchups_estimated_vs_observed_45sizeclasses_all.png'),'Resolution',600)
