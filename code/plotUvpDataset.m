% ======================================================================= %
%                                                                         %
% This script plots the POC flux data derived from the UVP5 dataset. It   %
% has 3 sections:                                                         %
%   Section 1 - Presets.                                                  %
%   Section 2 - Plot number of casts.                                     %
%   Section 3 - Plot Figure S3 (the UVP5-derived POC fluxes compared to   %
%               our trap and radionuclide compilation of POC fluxes).     %
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

%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Filename declarations 
filenameUvpProcessedDataset45sc  = 'pocflux_bisson_45sc_monthly_and_annual_all_depths.mat';
filenameMonthlyPocFlux           = 'pocflux_compilation_monthly.mat';
filenameTimeseriesInformation    = 'timeseries_station_information.mat';

% Define parameters
MOLAR_MASS_CARBON = 12.011; % g mol-1

% Load the UVP5 data
load(fullfile('.','data','processed','UVP5',filenameUvpProcessedDataset45sc),...
    'uvpFluxByCastAvg','uvpFluxByCastErr','uvpMonthlyFlux','targetDepths',...
    'castMonthlyDistrib')

maxNumCastsPerMonth = max(castMonthlyDistrib,[],'all');
nTargetDepths = numel(targetDepths);

% Load the trap and radionuclide compilation
load(fullfile('.','data','processed',filenameMonthlyPocFlux),...
    'obsRawProfileValues','obsRawProfileDataType','obsRawProfileDepths')

% Load station information
load(fullfile('.','data','processed',filenameTimeseriesInformation))

% Labels used for the plots
monthsLabel = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

%%
% -------------------------------------------------------------------------
% SECTION 2 - PLOT NUMBER OF CASTS
% -------------------------------------------------------------------------

categoricalLabelMonths = categorical(monthsLabel);
categoricalLabelMonths = reordercats(categoricalLabelMonths,monthsLabel);

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.42 0.27],'Color','w') 
y = castMonthlyDistrib(:,:);
h = bar(categoricalLabelMonths,y);

% Change the colors of each bar
coloursBars = jet(size(h,2)); 
coloursBars = mat2cell(coloursBars,ones(size(coloursBars,1),1),3);
set(h,{'FaceColor'},coloursBars)

box on
yl = ylabel('Number of UVP5 casts','FontSize',12);
lg = legend(STATION_NAMES);
set(lg,'Box','off') 

exportgraphics(gcf,fullfile('.','figures','uvp_numbercasts_45sizeclasses.png'),'Resolution',600)

%%
% -------------------------------------------------------------------------
% SECTION 3 - PLOT FIGURE S3 (THE UVP5-DERIVED POC FLUXES COMPARED TO OUR
% TRAP AND RADIONUCLIDE COMPILATION OF POC FLUXES)
% -------------------------------------------------------------------------

iSelectedDepths = 10:nTargetDepths;

for iLoc = 1:NUM_LOCS
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.25 0.75],'Color','w')
    haxis = zeros(12,1);

    isCastPlotted = 0;
    
    for iMonth = 1:12

        haxis(iMonth) = subaxis(4,3,iMonth,'Spacing',0.055,'Padding',0,'Margin',0.13);
        ax(iMonth).pos = get(haxis(iMonth),'Position');
        
        % Calculate xlim
        maxAnnualEst = max(uvpFluxByCastAvg(:,:,:,iLoc),[],'all','omitnan'); % mg C m-2 d-1
        maxAnnualObs = max(obsRawProfileValues(:,:,iLoc).*MOLAR_MASS_CARBON,[],'all','omitnan'); % mmol C m-2 d-1 --> mg C m-2 d-1 

        maxMonthlyObs = max(obsRawProfileValues(:,iMonth,iLoc).*MOLAR_MASS_CARBON,[],'omitnan');
        if (isnan(maxMonthlyObs)) 
            maxMonthlyObs = 0;
        end
        maxMonthlyEst = max(uvpFluxByCastAvg(:,:,iMonth,iLoc),[],'all','omitnan');
        if (isnan(maxMonthlyEst)) 
            maxMonthlyEst = 0;
        end

%         if (iLoc == 2 && iMonth ~= 8)
%             maxFlux = 190;
%         elseif (iLoc == 2 && iMonth == 8)
%             maxFlux = maxMonthlyEst;
%         else
            maxFlux = max(maxAnnualObs,maxAnnualEst);
%         end
        
        % Calculate ylim
        maxDepth = targetDepths(end);

        % Plot estimates from the UVP5: casts
        valse = squeeze(uvpFluxByCastAvg(:,iSelectedDepths,iMonth,iLoc));
        valse(valse==0) = NaN;
        valseresh = reshape(valse.', [], 1);
        ze = targetDepths(iSelectedDepths);
        zerep = repmat(ze,[maxNumCastsPerMonth 1]);
%         d1 = plot(haxis(iMonth),valseresh,zerep,'Color','r','LineWidth',0.4,...
%             'DisplayName','Estimated (individual casts)');
        % Replace red line by red scatters
        d1 = scatter(haxis(iMonth),valseresh,zerep,10,'r','filled',...
            'DisplayName','Estimated (individual casts)');
        
        hold on
                
        % Plot estimates from the UVP5: monthly mean
        valse = squeeze(uvpMonthlyFlux(iSelectedDepths,iMonth,iLoc,1));
        valse(valse==0) = NaN;
        ze = targetDepths(iSelectedDepths);
        d2 = plot(haxis(iMonth),valse,ze,'Color','k','LineWidth',1.5,...
            'DisplayName','Estimated (mean)');
        hold on
        
        % Plot observations
        
        % First - sediment trap
        mask = strcmp(squeeze(obsRawProfileDataType(:,iMonth,iLoc)),'trap');
        if (sum(mask) == 0)
            d3 = plot(NaN,NaN,'o','MarkerEdgeColor',[0.85 0.85 0.85],...
                'MarkerFaceColor',[0.85 0.85 0.85],'LineWidth',0.5,...
                'DisplayName','Observed (trap)');
        else
            d3 = plot(obsRawProfileValues(mask,iMonth,iLoc).*MOLAR_MASS_CARBON,...
                obsRawProfileDepths(mask,iMonth,iLoc),'o','MarkerEdgeColor',[0.85 0.85 0.85],...
                'MarkerFaceColor',[0.85 0.85 0.85],'LineWidth', 0.5,... 
                'DisplayName','Observed (trap)');
        end
        hold on
        
        % Second - radionuclides
        mask = strcmp(squeeze(obsRawProfileDataType(:,iMonth,iLoc)),'radionuclide');
        if (sum(mask) == 0)
            d4 = plot(NaN,NaN,'+','MarkerEdgeColor',[0.6 0.6 0.6],...
                'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',1.5,...
                'DisplayName', 'Observed (radionuclide)');
        else
            d4 = plot(obsRawProfileValues(mask,iMonth,iLoc).*MOLAR_MASS_CARBON,...
                obsRawProfileDepths(mask,iMonth,iLoc),'+','MarkerEdgeColor',[0.6 0.6 0.6],...
                'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',1.5,...
                'DisplayName','Observed (radionuclide)');
        end
        hold off
            
%         valso = squeeze(obsRawProfileValues(:,iMonth,iLoc)).*MOLAR_MASS_CARBON;
%         valso(valso==0) = NaN;
%         zo = squeeze(obsRawProfileDepths(:,iMonth,iLoc));
%         plot(haxis(iMonth),valso,zo,'o','MarkerEdgeColor',[0.85 0.85 0.85],'MarkerFaceColor',[0.8 0.8 0.8])
%         hold off
        
        box on
        
        if (iLoc == 1)
            monthlyMaxFlux = 500;
            xlim([0 monthlyMaxFlux])
            xticks([0 200 400])
            xticklabels({'0','200','400'})
        elseif (iLoc == 2) 
            monthlyMaxFlux = 200;
            xlim([0 monthlyMaxFlux])
            xticks([0 75 150])
            xticklabels({'0','75','150'})
        elseif (iLoc == 3 && iMonth ~= 5)
            monthlyMaxFlux = 250;
            xlim([0 monthlyMaxFlux])
            xticks([0 100 200])
            xticklabels({'0','100','200'})
        elseif (iLoc == 3 && iMonth == 5)
            monthlyMaxFlux = 800;
            xlim([0 monthlyMaxFlux])
            xticks([0 400 800])
            xticklabels({'0','400','800'})
        elseif (iLoc == 4 && iMonth == 2)
            monthlyMaxFlux = 400;
            xlim([0 monthlyMaxFlux])
            xticks([0 200 400])
            xticklabels({'0','200','400'})
        elseif (iLoc == 4 && iMonth == 12)
            monthlyMaxFlux = 1000;
            xlim([0 monthlyMaxFlux])
            xticks([0 400 800])
            xticklabels({'0','400','800'})
        elseif (iLoc == 4 && (iMonth ~= 12 || iMonth ~= 2))
            monthlyMaxFlux = 400;
            xlim([0 monthlyMaxFlux])
            xticks([0 200 400])
            xticklabels({'0','200','400'})
        elseif (iLoc == 5 && iMonth ~= 8)
            monthlyMaxFlux = 100;
            xlim([0 monthlyMaxFlux])
            xticks([0 50 100])
            xticklabels({'0','50','100'})
        elseif (iLoc == 5 && iMonth == 8)
            monthlyMaxFlux = 200;
            xlim([0 monthlyMaxFlux])
            xticks([0 100 200])
            xticklabels({'0','100','200'})
        elseif (iLoc == 6 && iMonth == 6)
            monthlyMaxFlux = 150;
            xlim([0 monthlyMaxFlux])
            xticks([0 70 140])
            xticklabels({'0','70','140'})
        elseif (iLoc == 6 && iMonth == 7)
            monthlyMaxFlux = 900;
            xlim([0 monthlyMaxFlux])
            xticks([0 400 800])
            xticklabels({'0','400','800'})
        elseif (iLoc == 6 && (iMonth ~= 6 || iMonth ~= 7))
            monthlyMaxFlux = 70;
            xlim([0 monthlyMaxFlux])
            xticks([0 25 50])
            xticklabels({'0','25','50'})
        end

        ylim([0 1000])
        yticks([0 250 500 750 1000])
        if (iMonth == 1 || iMonth == 4 || iMonth ==7 || iMonth == 10)
            yticklabels({'0','250','500','750','1000'})
        else
            yticklabels([])
        end
        axh = gca;
        axh.YAxis.TickDirection = 'out';
        axh.TickLength = [0.03, 0.03]; % make tick marks longer
        
        title(monthsLabel(iMonth),'FontSize',14)

        set(gca,'YDir','Reverse','XAxisLocation','Bottom','xlabel',[],'ylabel',[])
        
        if (castMonthlyDistrib(iMonth,iLoc) > 0 && isCastPlotted == 0 && iMonth > 2)
            lg = legend([d1 d2 d3 d4],'NumColumns',2);
            lg.Position(1) = 0.12; lg.Position(2) = 0.02;
            lg.ItemTokenSize = [15,1];
            lg.FontSize = 11;
            set(lg,'Box','off') 
            isCastPlotted = 1;
        end

    end % iMonth

    ax(1).pos(1) = ax(1).pos(1)+0.040; ax(1).pos(2) = ax(1).pos(2)+0.040; 
    ax(2).pos(1) = ax(2).pos(1)+0.025; ax(2).pos(2) = ax(2).pos(2)+0.040; 
    ax(3).pos(1) = ax(3).pos(1)+0.010; ax(3).pos(2) = ax(3).pos(2)+0.040;  
    ax(4).pos(1) = ax(4).pos(1)+0.040; ax(4).pos(2) = ax(4).pos(2)+0.040; 
    ax(5).pos(1) = ax(5).pos(1)+0.025; ax(5).pos(2) = ax(5).pos(2)+0.040; 
    ax(6).pos(1) = ax(6).pos(1)+0.010; ax(6).pos(2) = ax(6).pos(2)+0.040; 
    ax(7).pos(1) = ax(7).pos(1)+0.040; ax(7).pos(2) = ax(7).pos(2)+0.040; 
    ax(8).pos(1) = ax(8).pos(1)+0.025; ax(8).pos(2) = ax(8).pos(2)+0.040; 
    ax(9).pos(1) = ax(9).pos(1)+0.010; ax(9).pos(2) = ax(9).pos(2)+0.040; 
    ax(10).pos(1) = ax(10).pos(1)+0.040; ax(10).pos(2) = ax(10).pos(2)+0.040; 
    ax(11).pos(1) = ax(11).pos(1)+0.025; ax(11).pos(2) = ax(11).pos(2)+0.040; 
    ax(12).pos(1) = ax(12).pos(1)+0.010; ax(12).pos(2) = ax(12).pos(2)+0.040; 
    for iMonth = 1:12
        set(haxis(iMonth),'Position',ax(iMonth).pos) 
    end
    
    % Give common xlabel, ylabel and title to your figure
    % Create a new axis
    a = axes;
    t = title(STATION_NAMES(iLoc),'FontSize',16);
    xl = xlabel('POC flux (mg C m^{-2} d^{-1})','FontSize',16); % ,'FontWeight','bold'
    yl = ylabel('Depth (m)','FontSize',16);
    % Specify visibility of the current axis as 'off'
    a.Visible = 'off';
    % Specify visibility of Title, XLabel, and YLabel as 'on'
    t.Visible = 'on';
    xl.Visible = 'on';
    yl.Visible = 'on';
    yl.Position(1) = yl.Position(1) - 0.005; yl.Position(2) = yl.Position(2) + 0.02; 
    xl.Position(1) = 0.5; xl.Position(2) = xl.Position(2) + 0.065;
    t.Position(1) = t.Position(1); t.Position(2) = t.Position(2) + 0.035;

    exportgraphics(gcf,fullfile('.','figures',...
        strcat('monthly_flux_vs_depth_45sizeclasses_',STATION_TAGS{iLoc},'.png')),'Resolution',600)
  
end % iLoc 

