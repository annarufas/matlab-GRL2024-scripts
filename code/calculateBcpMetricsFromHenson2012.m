function [henson2012] = calculateBcpMetricsFromHenson2012(...
    isGlobalCalculation,listLocalLats,listLocalLons)

% CALCULATEBCPMETRICSFROMHENSON2012 Calculates five BCP metrics: Martin's b, 
% PEeff and Teff from 100 to 2000 m using the statistical fit of Henson et 
% al. (2012), as well as z* and Teff 100 to 1000 m as a function of b. The 
% five metrics are calculated using the annual means of the variables they 
% have dependencies on.
%
%   INPUT: 
%       isGlobalCalculation - choice (1 or 0)
%       listLocalLats       - list of local latitudes
%       listLocalLons       - list of local longitudes
%
%   OUTPUT:
%       isGlobalCalculation = 1, useful for visualising the global output 
%           of the Henson et al. (2012) algorithm. The output is a global 
%           array of annual mean values of the five BCP metrics, with no 
%           standard deviations.
%       isGlobalCalculation = 0 produces local annual mean values and 
%           standard deviations of the five BCP metrics.
% 
%   This script uses these external functions:
%       generateMCparameters.m - from FileExchange
%       propagateErrorWithMC.m - from FileExchange
%
%   ORIGINALLY WRITTEN BY S. HENSON, NOCS (2011)
%   MODIFIED BY A. RUFAS, UNIVERISTY OF OXFORD (2024)
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   This version - Completed 17 April 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Define parameters of the fitted PEeff-T curve

% The std is provided in the caption of Fig. 1 in Henson et al. (2011),
% PEeff = 0.23*exp(-0.08*SST)

fitPar.ave = [0.23, 0.08];
fitPar.stdev = [0.04, 0.01]; 
fitPar.stdevlow = [fitPar.ave(1)-fitPar.stdev(1), fitPar.ave(2)-fitPar.stdev(2)];
fitPar.stdevupp = [fitPar.ave(1)+fitPar.stdev(1), fitPar.ave(2)+fitPar.stdev(2)];

%% Define supplementary functions for Teff 100 to 1000 m and z*

% Notice that Henson et al. (2012) provide Teff from 100 to 200 m, so Teff 
% from 100 to 1000 m needs to be calculated, and we use Martin's b for
% that.

z0 = 100;
z1 = 1000;
funcTeff  = @(x) (z1/z0).^(-x); % x is Martin b
funcZstar = @(x) -(z1-z0)./log((z1/z0).^(-x)); % x is Martin b

%% Load data inputs for Henson's algorithm

[qSstMonthly,qChlaMonthly,qPar0Monthly] = extractHensonForcingDataForChosenLocations(...
    isGlobalCalculation,listLocalLats,listLocalLons);

%% Calculations

if isGlobalCalculation

    % Carr (2002) algorithm for NPP
    isZeuCarr = 0; % 1=Carr 2002, 0=Behrenfeld & Falkowski 1997 (as in Henson et al. 2012)
    qNppCarrMonthly = Carr2002algorithm(qChlaMonthly,qSstMonthly,qPar0Monthly,isZeuCarr); % mg C m-2 d-1

%     % Save the NPP array created
%     npp_avg = qNppCarrMonthly; 
%     npp_lon = (-179.5:1:179.5)';
%     npp_lat = (-89.5:1:89.5)';
%     save(fullfile('.','data','raw','npp_carr2002_global.mat'),'npp_avg','npp_lon','npp_lat')
    
    % Apply Henson's algorithm to calculate b, PEeff and Teff 100 to 2000 m
    % No error propagation, takes time
    qNppCarrAnnualMean = mean(qNppCarrMonthly(:,:,:),3,'omitnan');
    qNppCarrAnnualStd  = std(qNppCarrMonthly(:,:,:),0,3,'omitnan');
    qSstAnnualMean     = mean(qSstMonthly(:,:,:),3,'omitnan');
    
    [nr_chla,nc_chla,~] = size(qChlaMonthly);
    martinb  = NaN(nr_chla,nc_chla);
    teff2000 = NaN(size(martinb));
    peeff    = NaN(size(martinb));

    for iCol = 1:nc_chla
        for iRow = 1:nr_chla
            if (qNppCarrAnnualMean(iRow,iCol) > 0 && ~isnan(qNppCarrAnnualMean(iRow,iCol)))
                svi = qNppCarrAnnualStd(iRow,iCol)./qNppCarrAnnualMean(iRow,iCol);
            else
                svi = 0;
            end
            parArray = [qSstAnnualMean(iRow,iCol);fitPar.ave(1);fitPar.ave(2);qNppCarrAnnualMean(iRow,iCol);svi];
            peeff(iRow,iCol) = Henson2012peeff(parArray);  
            teff2000(iRow,iCol) = Henson2012teff2000(parArray);
            martinb(iRow,iCol) = Henson2012martinb(parArray);
        end
    end
    
    % Calculate z* and Teff 100 to 1000 m
    zstar = arrayfun(funcZstar,martinb);
    teff1000 = arrayfun(funcTeff,martinb);

    % Save arrays
    henson2012.martinb  = martinb;
    henson2012.peeff    = peeff;
    henson2012.teff2000 = teff2000;
    henson2012.teff1000 = teff1000;
    henson2012.zstar    = zstar;

    save(fullfile('.','data','processed','bcpmetrics_henson2012_global.mat'),'henson2012*')

else % local data
    
    nLocs = length(listLocalLats);
    trueMartinb   = NaN(nLocs,1);
    trueTeff2000  = NaN(size(trueMartinb));
    truePeeff     = NaN(size(trueMartinb));
    estimMartinb  = NaN(size(trueMartinb));
    estimTeff2000 = NaN(size(trueMartinb));
    estimPeeff    = NaN(size(trueMartinb));
    
     % Calculate local NPP from the Carr (2002) algorithm
    isZeuCarr = 0; % 1=Carr 2002, 0=Behrenfeld & Falkowski 1997 (as in Henson et al. 2012)
    qNppCarrMonthly = Carr2002algorithm(qChlaMonthly,qSstMonthly,qPar0Monthly,isZeuCarr); % mg C m-2 d-1

    for iLoc = 1:nLocs

        % Apply Hensons's algorithm to calculate b, PEeff and Teff 100 to 2000 m
        qSstAnnualMean = mean(qSstMonthly(iLoc,:,:),3,'omitnan');
        qSstAnnualStd = std(qSstMonthly(iLoc,:,:),0,3,'omitnan');
        qNppCarrAnnualMean = mean(qNppCarrMonthly(iLoc,:,:),3,'omitnan');
        qNppCarrAnnualStd = std(qNppCarrMonthly(iLoc,:,:),0,3,'omitnan');
        
        if (qNppCarrAnnualMean > 0 && ~isnan(qNppCarrAnnualMean))
            svi = qNppCarrAnnualStd/qNppCarrAnnualMean; % seasonal variation index of NPP
        else
            svi = 0;
        end
        
        parArray = [qSstAnnualMean;fitPar.ave(1);fitPar.ave(2);qNppCarrAnnualMean;svi];
        truePeeff(iLoc) = Henson2012peeff(parArray);  
        trueTeff2000(iLoc) = Henson2012teff2000(parArray);
        trueMartinb(iLoc) = Henson2012martinb(parArray);
            
        % Error propagation using MC
        A = generateMCparameters('gaussian',[qSstAnnualMean,qSstAnnualStd],'plot',false);
        B = generateMCparameters('gaussian',[fitPar.ave(1),fitPar.stdev(1)],'plot',false);
        C = generateMCparameters('gaussian',[fitPar.ave(2),fitPar.stdev(2)],'plot',false);
        D = generateMCparameters('gaussian',[qNppCarrAnnualMean,qNppCarrAnnualStd],'plot',false);
        D(D<0) = min(D(D>0));
        E = repmat(svi,[1 length(A)]);
        paramMatrix = [A;B;C;D;E];
        [midval_peeff,ci_peeff,funvals_peeff] = propagateErrorWithMC(@Henson2012peeff,paramMatrix,'plot',false);  
        [midval_teff2000,ci_teff2000,funvals_teff2000] = propagateErrorWithMC(@Henson2012teff2000,paramMatrix,'plot',false);
        [midval_martinb,ci_martinb,funvals_martinb] = propagateErrorWithMC(@Henson2012martinb,paramMatrix,'plot',false);
    
%         if (any(isnan(funvals_martinb)) || any(isinf(funvals_martinb)) || any(imag(funvals_martinb) ~= 0))
%             print('funvals contains anomalous values')
%         end

        estimMartinb(iLoc)  = midval_martinb;
        estimPeeff(iLoc)    = midval_peeff;
        estimTeff2000(iLoc) = midval_teff2000;

        % Calculate z* and Teff 100 to 1000 m as a function of b 
        if (~isnan(midval_martinb))
            
%             % The following three methods produce similar results. I am
%             % going to use Method C as it requires less steps
%             
%             % Method A
%             A = generateMCparameters('bootstrapDistribution',funvals_martinb);
%             [midval_teff1000_1,ci_teff1000_1,funvals_teff1000_1] = propagateErrorWithMC(funcTeff,A,'plot',false);
% 
%             % Method B - caveat: the error distribution is not gaussian
%             B = generateMCparameters('gaussian',[midval_martinb,(ci_martinb(2)-ci_martinb(1))/2]);
%             [midval_teff1000_2,ci_teff1000_2,funvals_teff1000_2] = propagateErrorWithMC(funcTeff,B,'plot',false);
%             
%             % Method C
%             [midval_teff1000_3,ci_teff1000_3,funvals_teff1000_3] = propagateErrorWithMC(funcTeff,funvals_martinb,'plot',false);
            
            [midval_teff1000,ci_teff1000,funvals_teff1000] = propagateErrorWithMC(funcTeff,funvals_martinb,'plot',false);
            [midval_zstar,ci_zstar,funvals_zstar] = propagateErrorWithMC(funcZstar,funvals_martinb,'plot',false);
        end

        if (nLocs < 10)
            fprintf('\nLocation %d',iLoc)
            fprintf('\nThe true b is %4.3f and the estimated b is %4.3f',trueMartinb(iLoc),estimMartinb(iLoc))
            fprintf('\nThe bounds estimated for b are %4.3f to %4.3f',ci_martinb(1),ci_martinb(2))
            fprintf('\nThe true PEeff is %4.3f and the estimated PEeff is %4.3f',truePeeff(iLoc),estimPeeff(iLoc))
            fprintf('\nThe bounds estimated for PEeff are %4.3f to %4.3f',ci_peeff(1),ci_peeff(2))
%             fprintf('\nThe estimated Teff 100 to 1000 m Method A is %4.3f, with bounds %4.3f to %4.3f',midval_teff1000_1,ci_teff1000_1(1),ci_teff1000_1(2))
%             fprintf('\nThe estimated Teff 100 to 1000 m Method B is %4.3f, with bounds %4.3f to %4.3f',midval_teff1000_2,ci_teff1000_2(1),ci_teff1000_2(2))
%             fprintf('\nThe estimated Teff 100 to 1000 m Method C is %4.3f, with bounds %4.3f to %4.3f\n',midval_teff1000_3,ci_teff1000_3(1),ci_teff1000_3(2))
            fprintf('\nThe estimated Teff 100 to 1000 m is %4.3f, with bounds %4.3f to %4.3f',midval_teff1000,ci_teff1000(1),ci_teff1000(2))
            fprintf('\nThe estimated z* is %4.0f, with bounds %4.0f to %4.0f\n',midval_zstar,ci_zstar(1),ci_zstar(2))
        end
        
        % Save arrays
        henson2012.martinb.ave(iLoc)       = midval_martinb;
        henson2012.martinb.stdevupp(iLoc)  = ci_martinb(2);
        henson2012.martinb.stdevlow(iLoc)  = ci_martinb(1);
        henson2012.martinb.max(iLoc)       = max(funvals_martinb);
        henson2012.martinb.min(iLoc)       = min(funvals_martinb);

        henson2012.teff1000.ave(iLoc)      = midval_teff1000;
        henson2012.teff1000.stdevupp(iLoc) = ci_teff1000(2);
        henson2012.teff1000.stdevlow(iLoc) = ci_teff1000(1);
        henson2012.teff1000.max(iLoc)      = max(funvals_teff1000);
        henson2012.teff1000.min(iLoc)      = min(funvals_teff1000);

        henson2012.zstar.ave(iLoc)         = midval_zstar;
        henson2012.zstar.stdevupp(iLoc)    = ci_zstar(2);
        henson2012.zstar.stdevlow(iLoc)    = ci_zstar(1);
        henson2012.zstar.max(iLoc)         = max(funvals_zstar);
        henson2012.zstar.min(iLoc)         = min(funvals_zstar);

    end % iLoc
    
    save(fullfile('.','data','processed','bcpmetrics_henson2012_local.mat'),'henson2012*')

end % isGlobalCalculation

% %% Contour plot of the variables controlling b
% 
% xtemp         = zeros(nr_ssta*nc_ssta,1); % for SST
% ynpp          = zeros(size(xtemp)); % for NPP
% zb            = zeros(size(xtemp)); % for b
% varflux2000   = zeros(size(xtemp));
% varexportprod = zeros(size(xtemp));
% 
% idx = 1;
% for iCol = 1:nc_ssta
%     for iRow = 1:nr_ssta
%         xtemp(idx) = qSstAnnualMean(iRow,iCol);
%         ynpp(idx) = qNppCarrAnnualMean(iRow,iCol);
%         zb(idx) = henson2012martinb(iRow,iCol);
%         varflux2000(idx) = henson2012ep2000(iRow,iCol);
%         varexportprod(idx) = henson2012ep100(iRow,iCol);
%         idx = idx + 1;
%     end
% end
% 
% [xtemp_sort, sortIdx] = sort(xtemp);
% ynpp_sort = ynpp(sortIdx);
% zb_sort = zb(sortIdx);
% 
% idxFirstNanTemp = find(isnan(xtemp_sort), 1, 'first');
% crop_xtemp_sort = xtemp_sort(1:idxFirstNanTemp-1);
% crop_ynpp_sort = ynpp_sort(1:idxFirstNanTemp-1);
% crop_ynpp_sort(isnan(crop_ynpp_sort)) = 0;
% crop_zb_sort = zb_sort(1:idxFirstNanTemp-1);
% crop_zb_sort(isnan(crop_zb_sort)) = 0;
% 
% mintemp = min(crop_xtemp_sort(:));
% maxtemp = max(crop_xtemp_sort(:));
% minnpp = min(crop_ynpp_sort(:));
% maxnpp = max(crop_ynpp_sort(:));
% 
% xv = linspace(mintemp, maxtemp, 150);
% yv = linspace(minnpp, maxnpp, 150);
% [Xm,Ym] = ndgrid(xv, yv);
% Zm = griddata(crop_xtemp_sort, crop_ynpp_sort, crop_zb_sort, Xm, Ym);
% 
% figure()
% set(gcf,'Units','Normalized','Position',[0.01 0.05 0.25 0.35],'Color','w') 
% contourf(Xm,Ym,Zm) % ,'ShowText','on'
% ylim([0 1000])
% xlim([-2.5 29])
% xh = xlabel('SST (ºC)');
% yh = ylabel('NPP (mg C m^{-2} d^{-1})');
% th = title('Henson et al. (2012) Martin b','FontSize',12);
% th.Position(2) = th.Position(2) + 30;
% xh.Position(2) = xh.Position(2) - 20;            
% yh.Position(1) = yh.Position(1) - 1;
% colorbar
% grid
% set(haxis(iSubplot),'FontSize',12)
% set(gcf,'PaperPositionMode','auto')
% print(gcf,fullfile(fullpathPlotsDir,'henson2012countour'),'-dpdf','-r0')  
% iFig = iFig + 1;
% 
% % [zb_sort2, sortIdx2] = sort(crop_zb_sort);
% % xtemp_sort2          = crop_xtemp_sort(sortIdx2);
% % zb_sort2             = crop_zb_sort(sortIdx2);
% % 
% % idxLastZeroZ       = find(zb_sort2 == 0, 1, 'last');
% % crop_xtemp_sort2   = xtemp_sort2(idxLastZeroZ+1:end);
% % crop_zb_sort2      = zb_sort2(idxLastZeroZ+1:end);
% % 
% % mdl = fitlm(crop_xtemp_sort2,crop_zb_sort2);
% % 
% % mymat = exportProduction(:,:);
% % integratedPOC = integrateAnnualPocAccordingToAreaOfTheGridCell(WOA13_lat,WOA13_lon,mymat);
% 
% end
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

% *************************************************************************

function [peeff] = Henson2012peeff(x)
    
    param.sst = x(1);
    param.sst_coeff1 = x(2);
    param.sst_coeff2 = x(3);

    % Particle export efficiency (Eq. S1 in the Supplementary Material of Henson et al. 2012)
    if (~isnan(param.sst))
        peeff = param.sst_coeff1*exp(-param.sst_coeff2*param.sst); % 0-1
    else
        peeff = 0;
    end
    
end % Henson2012peeff

% *************************************************************************

function [teff2000] = Henson2012teff2000(x)

    param.sst = x(1);
    param.sst_coeff1 = x(2);
    param.sst_coeff2 = x(3);
    param.npp = x(4);
    param.svi = x(5);

    % Export production
    if (~isnan(param.npp))
        ep = param.npp*Henson2012peeff(x); % mg m-2 d-1
    else
        ep = 0;
    end

    % The calculation of POC flux at 2000 m needs prd, rld and prr
    if (~isnan(param.svi))
        prd = 1e-3*(31*param.svi^2 + 49*param.svi + 7.8);
        rld = 1400*exp(-0.54*param.svi);
        prr = 1e-3*(2.6*param.svi^2 - 4.2*param.svi + 4.8);
    else
        prd = 0;
        rld = 0;
        prr = 0;
    end
    if (~isnan(param.npp))
        pocFlux2000 = param.npp*((prd*exp(-(2000-100)/rld))+prr); % mg m-2 d-1
    else
        pocFlux2000 = 0;
    end

    % Teff surface --> 2000 m
    if (ep > 0)
        teff2000 = pocFlux2000/ep;
    else
        teff2000 = 0;
    end

end % Henson2012teff2000

% *************************************************************************

function [martinb] = Henson2012martinb(x)

    teff2000 = Henson2012teff2000(x);

    % Martin's b, based on the power-law function F_2000 = EP x (2000/100)^b
    if (teff2000 > 0 && ~isnan(teff2000))
        martinb = -1.*(log(teff2000)/log(2000/100));
    else
        martinb = 0;
    end
    
end % Henson2012martinb

% *************************************************************************

function [npp] = Carr2002algorithm(chla, sst, par0, isZeuCarr)
    
    % Calculate light-related properties

    [nr_chla,nc_chla,nt_chla] = size(chla);
    kPar = NaN(nr_chla,nc_chla,nt_chla); % attenuation coefficient for PAR 
    parZeu = NaN(nr_chla,nc_chla,nt_chla); % PAR at the euphotic layer depth 
    zeuCarr = NaN(nr_chla,nc_chla,nt_chla); % euphotic layer depth, Carr (2002) 
    zeuBF = NaN(nr_chla,nc_chla,nt_chla); % euphotic layer depth,  Behrenfeld & Falkowski (1997) 

    for iRow = 1:nr_chla
        for iCol = 1:nc_chla
            for iMonth = 1:12
                if (~isnan(chla(iRow,iCol,iMonth)))

                    % Carr (2002) Eq. 1, after Nelson and Smith (1991)
                    kPar(iRow,iCol,iMonth) = 0.04 + (0.0088.*chla(iRow,iCol,iMonth))... 
                        + 0.054.*chla(iRow,iCol,iMonth).^(0.66); % m-1 

                    % Carr (2002) Eq. 4
                    zeuCarr(iRow,iCol,iMonth) = -log(0.01)./kPar(iRow,iCol,iMonth); % m

                    % Henson et al. (2012) use the Behrenfeld & Falkowski (1997) 
                    % algorithm (VGPM model, https://sites.science.oregonstate.edu/ocean.productivity/vgpm.code.php)
                    % to compute zeu instead of Carr (2002) proposed equation
                    if (chla(iRow,iCol,iMonth) < 1)
                        Ctot=38*(chla(iRow,iCol,iMonth).^0.425); % integrated water column chl
                    else
                        Ctot=40.2*(chla(iRow,iCol,iMonth).^0.507);
                    end
                    zeu = 568.2*(Ctot.^-0.746);
                    if (zeu > 102)
                        zeu = 200*(Ctot.^-0.293);
                    end
                    zeuBF(iRow,iCol,iMonth) = zeu; % m

                    if (isZeuCarr)
                        zeu = zeuCarr;
                    else
                        zeu = zeuBF;
                    end

                    %Carr (2002) Eq. 2, after Riley (1957)
                    parZeu(iRow,iCol,iMonth) = par0(iRow,iCol,iMonth)...
                        .*(1-exp(-kPar(iRow,iCol,iMonth).*zeu(iRow,iCol,iMonth)))...
                        ./(kPar(iRow,iCol,iMonth).*zeu(iRow,iCol,iMonth)); % W m-2 

                end
            end
        end
    end

%     % Visual comparison of zeu from the Carr (2002) vs the B&F algorithm
%     
%     nSubplots = 2;
%     mytitlestring = {'Carr (2002) algorithm','B&F algorithm'};
%     figure()
%     set(gcf,'Units','Normalized','Position',[0.01 0.05 0.29 0.49],'Color','w')
%     haxis = zeros(nSubplots,1);
%     
%     for iSubplot = 1:nSubplots  
%         haxis(iSubplot) = subplot(2,1,iSubplot);
%         ax(iSubplot).pos = get(haxis(iSubplot),'Position');
%         if (iSubplot == 1)
%             mydata = zeuCarr;
%         elseif (iSubplot == 2)
%             mydata = zeuBF;
%         end
%         pcolor(flipud(rot90(mydata(:,:,4))))
%         caxis([40 110])
%         shading interp; colormap(jet)
%         set(gca,'xticklabels',[],'yticklabels',[])
%         box on
%         title(sprintf('%s',mytitlestring{iSubplot})) 
%         set(haxis(iSubplot),'FontSize',12)
%     end
%     ax(1).pos(1) = ax(1).pos(1) - 0.12; ax(1).pos(4) = ax(1).pos(4) + 0.040;
%     ax(2).pos(1) = ax(2).pos(1) - 0.12; ax(2).pos(4) = ax(2).pos(4) + 0.040; 
%     for iSubplot = 1:nSubplots
%         set(haxis(iSubplot), 'Position', ax(iSubplot).pos) 
%     end    
%     cb = colorbar('Location','eastoutside');
%     cb.Position(1) = 0.85;
%     cb.Position(2) = 0.25;
%     cb.Position(4) = 0.60; 
%     cb.Label.String = 'z_{eu} (m)'; 
%     cb.FontSize = 12;

    % .....................................................................
    
    % Calculate NPP (Eq. 3 in Carr 2002)

    % Two codes presented: as originally coded by me and as coded by S.
    % Henson. Both codes are equivalent, the difference is in the way they 
    % present the Epplye curve (max. phytoplankton growth rate as a 
    % function of temperature). Eppley (1972) originally formulated his 
    % equation for max. growth rate in units of d-1 (umax). Later on, in 
    % 1991, Platt et al. revisited Eppley (1972) equation so that it could 
    % be expressed in units of mg C (mg chla)-1 d-1 (Pmax). Platt et al. 
    % (1991) equation implicitly assumes a C/chla ratio of (=24/0.59) 
    % 40.1 mg C (mg chl)-1. In my calculations, I had originally assumed a 
    % C/chla of 50 mg C (mg chl)-1 (after Fasham 1990). 

%     % As originally coded by me (umax)
%     THETA_C_TO_CHL = 40.1; % 50 g C (g chla)-1 = 0.020 g chla (g C)-1 (Fasham et al. 1990)
%     muMax = 0.59.*exp(0.0633.*sst(:,:,:)); % d-1, the Eppley curve (Eppley 1972), after Kremer et al. (2017) – 0.59 d-1 is the max. growth rate at 0ºC
%     % muMax = 0.59.*1.88.^(sst(:,:,:)./10) % d-1, the Eppley curve (Eppley 1972), after Kremer et al. (2017) - Q10 factor equivalent of the exponential model
%     Ik = muMax(:,:,:)./(2.64/THETA_C_TO_CHL); % W m-2
%     chlaProd = chla(:,:,:).*((muMax(:,:,:).*parZeu(:,:,:))./(Ik(:,:,:)+parZeu(:,:,:))); % mg chla m-3 d-1
%     npp = chlaProd(:,:,:).*THETA_C_TO_CHL.*zeu(:,:,:); % mg chla m-3 d-1 --> mg C m-2 d-1

    % As coded by Steph (Pmax)
    ALPHA = 2.64; % mg C (mg chla)-1 d-1 (W m-2)-1
    Pmax = 24.*exp(0.09.*sst(:,:,:)); % mg C (mg chla)-1 d-1, the Eppley curve (Eppley 1972), after Platt et al. (1991) (see here: https://www.sciencedirect.com/science/article/pii/S0964274997800182?via=ihub)
    npp = chla(:,:,:).*(((Pmax(:,:,:).*parZeu(:,:,:))./((Pmax(:,:,:)./ALPHA)+parZeu(:,:,:)))).*zeu(:,:,:); % mg C m-2 d-1

%     % Compare with Fig. 4a in Henson et al. (2012)
%     nppAnnualMean(:,:) = (365/1000).*mean(npp(:,:,:),3,'omitnan'); % g C m-2 yr-1
%     figure(); pcolor(flipud(rot90(nppAnnualMean))); 
%     caxis([0 500]); 
%     cb = colorbar('FontSize', 15, 'FontWeight', 'bold'); 
%     cb.Label.String = 'NPP (g C m^{-2} yr^{-1})';
%     shading interp; colormap(jet); box on;

end % Carr2002algorithm

end