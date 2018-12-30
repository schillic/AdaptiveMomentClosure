classdef Constants
%CONSTANTS Constant definitions
    
    properties (Constant)
        % default value for event which does not change a parameter value
        C_n_defaultNoEvent = -1;
        
        % parameter domain types
        % standard domain
        C_si_domainStandard = '';
        % constant parameter whose value does not change during search
        C_si_domainConstant = 'c';
        % Boolean domain
        C_si_domainBoolean = 'b';
        
        % figure positions
        C_n_figureStandard = 0;
        C_n_figureZero = 1;
        C_n_figureGrid = 2;
        
        % mean/variance plot modes (must be a value x such that 0 <= x < 1)
        C_n_plotData = 0;
        C_n_plotSimulation = 0.1;
        C_n_plotBestClosure = 0.2;
        
        % plot modes
        C_si_plotLine = '';
        C_si_plotDot = '.';
        C_si_plotCircle = 'o';
        C_si_plotCross = 'x';
        C_si_plotPlus = '+';
        
        % plot window size
        C_n_figureBorderHeight = 27;
        C_n_plotWidth = 500;
        C_n_plotHeight = 500;
        C_n_plotParametersWidth = 500;
        C_n_plotParametersHeight = 500;
        
        % ignored measurements for ODE solver
        C_n_ignoreNone = 0;
        C_n_ignoreAll = 1;
        C_n_ignoreFirst = 2;
        
        % mode to compute the likelihood
        C_n_likelihoodSingle = 0;
        C_n_likelihoodTwo = 1;
    end
end