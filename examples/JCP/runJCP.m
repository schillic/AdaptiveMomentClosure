function runJCP( )
%% RUNJCP Run script for the JCP example

% reference parameters (for reference simulation)
c_paramsReference = cell(4, 1);
c_paramsReference{1} = QuickParameter(2*10^(-3), []); % a
c_paramsReference{2} = QuickParameter(2*10^(-1), []); % b
c_paramsReference{3} = QuickParameter(2*10^(-4), []); % c
c_paramsReference{4} = QuickParameter(2*10^(-3), []); % d
v_eventsReference = [];
o_paramsReference = ParameterVector(c_paramsReference, v_eventsReference);

% initial parameter estimation
c_paramsInit = cell(4, 1);
c_paramsInit{1} = QuickParameter(2*10^(-3) * 10, []); % a
c_paramsInit{2} = QuickParameter(2*10^(-1) * 10, []); % b
c_paramsInit{3} = QuickParameter(2*10^(-4) * 10, []); % c
c_paramsInit{4} = QuickParameter(2*10^(-3) * 2, []); % d
v_eventsInit = [];
o_paramsInit = ParameterVector(c_paramsInit, v_eventsInit);

% parameter domains
c_domains = {
    [], [], Constants.C_si_domainStandard; % a
    [], [], Constants.C_si_domainStandard; % b
    [], [], Constants.C_si_domainStandard; % c
    [], [], Constants.C_si_domainStandard; % d
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
    'ld', 4; ... % crashes
    };

% Boolean flags
s_flags = struct( ...
    );

% numbers
s_numbers = struct( ...
    'simulationItData', 5000, ... % number of simulations for reference data
    'simulationItLoop', 200, ... % number of simulations during search
    'iterations', 4000, ... % number of iterations / generated parameter samples
    'epsilon', 0.2, ... % size of the epsilon ball
    'newParameterVariance', 0.02, ... % variance of selection distribution
    'figurePlacement', Constants.C_n_figureGrid ... % figure placement
    );

% strings
s_strings = struct( ...
    'netFileName', 'JCP' ... % *.net file name (without file extension)
    );

s_vectors = struct( ...
    'measurementTimes', (50 : 50 : 500) ... % measurement time points
    );

s_cells = struct( ...
    'domains', {c_domains}, ... % parameter domains
    'measuredVariableNames', {{'P'}}, ... % measured species
    'closureMethods', {c_closureMethodsAvailable} ... % closure methods used
    );

s_objects = struct( ...
    'parametersReference', o_paramsReference, ... % true parameters
    'parametersInitial', o_paramsInit ... % initial parameter estimation
    );

% call general wrapper function with specified arguments
mainWrapper(s_flags, s_numbers, s_strings, s_vectors, s_cells, s_objects);