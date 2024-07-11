function [fitMartinCurve,fitExpCurve] = defineFittypesForPocFluxAttenuationCurves()

% DEFINEFITTYPESFORPOCFLUXATTENUATIONCURVES Define the functions containing 
% the coefficients b and z* and use Matlab's curve fitting toolbox along
% with a non-linear least squares method to determine the values of b and
% z*.
%
%   OUTPUT:
%       fitMartinCurve - resolves Martin's b
%       fitExpCurve    - resolves z*
%                              
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 9 April 2024   
%
% =========================================================================
%%

% I have chosen not to linearise the power-law and exponential curves
% (i.e., by applying log10) to treat every data point as equally important.

funcMartinCurve = @(b,f0,z0,z)     f0.*((z/z0).^(-b));      % Fz = F100 x (z/100)^b
funcExpCurve    = @(zstar,f0,z0,z) f0.*exp(-(z-z0)./zstar); % Fz = F100 x e^(-(z-100)/z*)
opts            = fitoptions('Method','NonlinearLeastSquares'); 

fitMartinCurve = fittype(funcMartinCurve,'problem',{'f0','z0'},...
    'independent','z','coefficients','b','options',opts);
fitExpCurve = fittype(funcExpCurve,'problem',{'f0','z0'},...
    'independent','z','coefficients','zstar','options',opts);

end