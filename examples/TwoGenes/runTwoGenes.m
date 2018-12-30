function runTwoGenes( )
%% RUNTWOGENES Run script for the TwoGenes example

% true parameters (for reference simulation)
c_paramsReference = cell(7, 1);
c_paramsReference{1} = QuickParameter(0.1, []); % c1
c_paramsReference{2} = QuickParameter(0.01, []); % c2M
c_paramsReference{3} = QuickParameter(10, []); % c3
c_paramsReference{4} = QuickParameter(0.005, []); % c4
c_paramsReference{5} = QuickParameter(0.1, []); % c5
c_paramsReference{6} = QuickParameter(10, []); % c6
c_paramsReference{7} = QuickParameter(0.1, []); % d
v_eventsReference = [];
o_paramsReference = ParameterVector(c_paramsReference, v_eventsReference);

% initial parameter estimation
c_paramsInit = cell(7, 1);
c_paramsInit{1} = QuickParameter(0.1 * 0.1, []); % c1
c_paramsInit{2} = QuickParameter(0.01 * 5, []); % c2M
c_paramsInit{3} = QuickParameter(10 * 0.1, []); % c3
c_paramsInit{4} = QuickParameter(0.005 * 5, []); % c4
c_paramsInit{5} = QuickParameter(0.1 * 0.1, []); % c5
c_paramsInit{6} = QuickParameter(10 * 0.1, []); % c6
c_paramsInit{7} = QuickParameter(0.1 * 0.1, []); % d
v_eventsInit = [];
o_paramsInit = ParameterVector(c_paramsInit, v_eventsInit);

% parameter domains
c_domains = {
    [], [], Constants.C_si_domainStandard; % c1
    [], [], Constants.C_si_domainStandard; % c2
    [], [], Constants.C_si_domainStandard; % c3
    [], [], Constants.C_si_domainStandard; % c4
    [], [], Constants.C_si_domainStandard; % c5
    [], [], Constants.C_si_domainStandard; % c6
    [], [], Constants.C_si_domainStandard % d
    };

% approximations (closures) available
c_closureMethodsAvailable = {...
    'dm', 2; ...
    'zc', 2; ... % crashes
    'zv', 2; ... % crashes
    'dm', 3; ...
    'zc', 3; ... % crashes
    'zv', 3; ... % crashes
    'ld', 3; ...
    'dm', 4; ...
    'zc', 4; ... % crashes
    'zv', 4; ... % crashes
    'ld', 4; ... % crashes
    };

% Boolean flags
s_flags = struct( ...
    );

% numbers
s_numbers = struct( ...
    'simulationItData', 5000, ... % number of simulations for reference data
    'simulationItLoop', 100, ... % number of simulations during search
    'iterations', 10000, ... % number of iterations / generated parameter samples
    'epsilon', 0.2, ... % size of the epsilon ball
    'newParameterVariance', 0.03, ... % variance of selection distribution
    'likelihoodMode', Constants.C_n_likelihoodTwo, ... % likelihood mode
    'figurePlacement', Constants.C_n_figureGrid ... % figure placement
    );

% strings
s_strings = struct( ...
    'netFileName', 'TwoGenes' ... % *.net file name (without file extension)
    );

s_vectors = struct( ...
    'measurementTimes', (10 : 10 : 100) ... % measurement time points
    );

s_cells = struct( ...
    'domains', {c_domains}, ... % parameter domains
    'measuredVariableNames', {{'P1', 'P2'}}, ... % measured species
    'closureMethods', {c_closureMethodsAvailable} ... % closure methods used
    );

s_objects = struct( ...
    'parametersReference', o_paramsReference, ... % true parameters
    'parametersInitial', o_paramsInit ... % initial parameter estimation
    );

% call general wrapper function with specified arguments
mainWrapper(s_flags, s_numbers, s_strings, s_vectors, s_cells, s_objects);