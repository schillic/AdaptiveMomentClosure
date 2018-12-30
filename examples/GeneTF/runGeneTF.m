function runGeneTF( )
%% RUNGENETF Run script for the GeneTF example

% true parameters (for reference simulation)
c_paramsReference = cell(6, 1);
c_paramsReference{1} = QuickParameter(2*10^(-2), []); % a
c_paramsReference{2} = QuickParameter(2*10^(-3), []); % b
c_paramsReference{3} = QuickParameter(2*10^(-4), []); % c
c_paramsReference{4} = QuickParameter(2*10^(-1), []); % d
c_paramsReference{5} = QuickParameter(2*10^(-4), []); % e
c_paramsReference{6} = QuickParameter(2*10^(-2), []); % f
v_eventsReference = [];
o_paramsReference = ParameterVector(c_paramsReference, v_eventsReference);

% initial parameter estimation
c_paramsInit = cell(6, 1);
c_paramsInit{1} = QuickParameter(2*10^(-2) * 10, []); % a
c_paramsInit{2} = QuickParameter(2*10^(-3) * 10, []); % b
c_paramsInit{3} = QuickParameter(2*10^(-4) * 10, []); % c
c_paramsInit{4} = QuickParameter(2*10^(-1) / 10, []); % d
c_paramsInit{5} = QuickParameter(2*10^(-4) * 10, []); % e
c_paramsInit{6} = QuickParameter(2*10^(-2) / 10, []); % f
v_eventsInit = [];
o_paramsInit = ParameterVector(c_paramsInit, v_eventsInit);

% parameter domains
c_domains = {
    [], [], Constants.C_si_domainStandard; % a
    [], [], Constants.C_si_domainStandard; % b
    [], [], Constants.C_si_domainStandard; % c
    [], [], Constants.C_si_domainStandard; % d
    [], [], Constants.C_si_domainStandard; % e
    [], [], Constants.C_si_domainStandard % f
    };

% approximations (closures) available
c_closureMethodsAvailable = {...
    'dm', 2; ...
    'zc', 2; ...
    'zv', 2; ... % crashes
    'dm', 3; ...
    'zc', 3; ...
    'zv', 3; ... % crashes
    'ld', 3; ...
    'dm', 4; ...
    'zc', 4; ...
    'zv', 4; ... % crashes
    'ld', 4; ...
    };

% Boolean flags
s_flags = struct( ...
    'computeDataFile', 0, ... % compute reference data?
    'computeClosureFiles', 0, ... % compute closure files?
    'computeSimulationFile', 0, ... % compute simulation script?
    'removeCrashingClosures', 0 ... % remove closures after crashing?
    );

% numbers
s_numbers = struct( ...
    'simulationItData', 5000, ... % number of simulations for reference data
    'simulationItLoop', 100, ... % number of simulations during search
    'iterations', 5000, ... % number of iterations / generated parameter samples
    'epsilon', 0.2, ... % size of the epsilon ball
    'newParameterVariance', 0.03, ... % variance of selection distribution
    'seed', 0, ... % random number generation seed (0 = default, < 0 = none)
    'rejectionLimit', 100, ... % reevaluate closure after a number of rejections
    'figurePlacement', Constants.C_n_figureGrid ... % figure placement
    );

% strings
s_strings = struct( ...
    'netFileName', 'GeneTF', ... % *.net file name (without file extension)
    'simulationName', 'GillespieGeneTF', ... % simulation function name
    'figureName', 'GeneTF' ... % figure name
    );

s_vectors = struct( ...
    'measurementTimes', (100 : 100 : 1000), ... % measurement time points
    'plotGrid', [3, 2] ... % plot grid
    );

s_cells = struct( ...
    'domains', {c_domains}, ... % parameter domains
    'measuredVariableNames', {{'P'}}, ... % measured species
    'funNamesOdeSolve', {{'sundials', {'matlab', 5}}}, ... % ODE solver function
    'closureMethods', {c_closureMethodsAvailable} ... % closure methods used
    );

s_objects = struct( ...
    'parametersReference', o_paramsReference, ... % true parameters
    'parametersInitial', o_paramsInit ... % initial parameter estimation
    );

% call general wrapper function with specified arguments
mainWrapper(s_flags, s_numbers, s_strings, s_vectors, s_cells, s_objects);