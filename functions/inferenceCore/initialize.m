function [ s_constants, s_protocol, s_current, s_statistics, s_closures, ...
    s_plotSetup ] = initialize( ...
        s_flags, s_numbers, s_strings, s_vectors, s_cells, s_objects )
%INITIALIZE Initialize the data structures

% start time measuring
timer_refTotal = tic;

% initialize internal data structure (with future constant terms) from input
s_constants = setupConstants(s_flags, s_numbers, s_strings, s_vectors, ...
    s_cells, s_objects);

% reset RNG
resetRng(s_constants.seed);

% initialize folder and file names
s_constants = setupFiles(s_constants);

% read the model
s_net = readNet(s_constants.netFile);

% initialize main remaining data structures
[s_protocol, s_current, s_constants] = ...
    setupDataStructures(s_constants, s_net.parameter);

% initialize parameter domains
s_constants = setupDomains(s_constants);

% initialize indices of the measured species (for reading simulation/ODE data)
s_constants = setupSpeciesIndices(s_constants, s_net.species);

% initialize ODE solving functions
s_constants = setupOdeSolvingFunctions(s_constants);

% generate simulation file
s_constants = generateSimulationFile(s_constants, s_net);

% reset RNG
resetRng(s_constants.seed);

% get reference data
s_constants = setupReferenceMoments(s_constants);

% reset RNG
resetRng(s_constants.seed);

% precompute the equations for each closure method
[s_constants, s_closures] = setupClosures(s_constants, s_net);

% general plot structure
s_plotSetup = setupPlotStructure();

% dynamic switches
s_constants = setupDynamicSwitches(s_constants);

% initialize statistics data structure
s_statistics = setupStatistics(timer_refTotal);

end
%% --- end of main function ---


%% --- helper functions ---

function checkSufficientInputs(s_flags, s_numbers, s_strings, s_vectors, ...
    s_cells, s_objects) %#ok<INUSL>
%% checks that all relevant information is provided by the user

    assert((...
        isfield(s_numbers, 'simulationItData') && ...
        isfield(s_numbers, 'simulationItLoop') && ...
        isfield(s_numbers, 'iterations') && ...
        isfield(s_numbers, 'epsilon') && ...
        isfield(s_numbers, 'newParameterVariance')), ...
        'A relevant number is missing.');
    assert((...
        isfield(s_strings, 'netFileName')), ...
        'A relevant string is missing.');
    assert((...
        isfield(s_vectors, 'measurementTimes')), ...
        'A relevant vector is missing.');
    assert((...
        isfield(s_cells, 'domains') && ...
        isfield(s_cells, 'measuredVariableNames') && ...
        isfield(s_cells, 'closureMethods')), ...
        'A relevant cell is missing.');
    assert((...
        isfield(s_objects, 'parametersReference') && ...
        isfield(s_objects, 'parametersInitial')), ...
        'A relevant object is missing.');
end


function [ s_constants ] = setupConstants( s_flags, s_numbers, s_strings, ...
        s_vectors, s_cells, s_objects )
%% packs all constant terms from the input data structures

    % optional parameters with default settings
    
    if (isfield(s_flags, 'computeDataFile'))
        b_computeDataFile = s_flags.computeDataFile;
    else
        b_computeDataFile = false;
    end
    if (isfield(s_flags, 'computeClosureFiles'))
        b_computeClosureFiles = s_flags.computeClosureFiles;
    else
        b_computeClosureFiles = false;
    end
    if (isfield(s_flags, 'computeSimulationFile'))
        b_computeSimulationFile = s_flags.computeSimulationFile;
    else
        b_computeSimulationFile = false;
    end
    if (isfield(s_flags, 'removeCrashingClosures'))
        b_removeCrashingClosures = s_flags.removeCrashingClosures;
    else
        b_removeCrashingClosures = false;
    end
    if (isfield(s_flags, 'plot'))
        b_plot = s_flags.plot;
    else
        b_plot = true;
    end
    if (isfield(s_numbers, 'likelihoodMode'))
        n_likelihoodMode = s_numbers.likelihoodMode;
    else
        n_likelihoodMode = Constants.C_n_likelihoodSingle;
    end
    if (isfield(s_numbers, 'rejectionLimit'))
        n_rejectionLimit = s_numbers.rejectionLimit;
    else
        n_rejectionLimit = -1;
    end
    if (isfield(s_numbers, 'likelihoodThreshold'))
        n_likelihoodThreshold = s_numbers.likelihoodThreshold;
    else
        n_likelihoodThreshold = -Inf;
    end
    if (isfield(s_numbers, 'figurePlacement'))
        n_figurePlacement = s_numbers.figurePlacement;
    else
        n_figurePlacement = Constants.C_n_figureStandard;
    end
    if (isfield(s_numbers, 'seed'))
        n_seed = s_numbers.seed;
    else
        n_seed = 0;
    end
    if (isfield(s_strings, 'baseName'))
        si_baseName = s_strings.baseName;
    else
        si_baseName = s_strings.netFileName;
    end
    if (isfield(s_strings, 'appendName'))
        si_appendName = s_strings.appendName;
    else
        si_appendName = '';
    end
    if (isfield(s_strings, 'simulationName'))
        si_simulationName = s_strings.simulationName;
    else
        si_simulationName = ''; % set later
    end
    if (isfield(s_strings, 'figureName'))
        si_figureName = s_strings.figureName;
    else
        si_figureName = si_baseName;
    end
    if (isfield(s_strings, 'simulationFile'))
        si_simulationFile = s_strings.simulationFile;
    else
        si_simulationFile = ''; % set later
    end
    if (isfield(s_strings, 'diaryName'))
        si_diaryName = s_strings.diaryName;
    else
        si_diaryName = [si_baseName, '-diary.txt'];
    end
    if (isfield(s_vectors, 'dynamicSwitches'))
        v_dynamicSwitches = s_vectors.dynamicSwitches;
    else
        v_dynamicSwitches = [];
    end
    if (isfield(s_vectors, 'momentTimes'))
        v_momentTimes = s_vectors.momentTimes;
    else
        v_momentTimes = []; % set later
    end
    if (isfield(s_cells, 'funNamesOdeSolve'))
        c_funNamesOdeSolve = s_cells.funNamesOdeSolve;
    else
        c_funNamesOdeSolve = {'sundials', {'matlab', 5}};
    end
    if (isfield(s_cells, 'timeFunctions'))
        c_timeFunctions = s_cells.timeFunctions;
    else
        c_timeFunctions = {};
    end
    
    % check that all mandatory inputs are present
    checkSufficientInputs(...
        s_flags, s_numbers, s_strings, s_vectors, s_cells, s_objects);
    
    
    s_constants = struct(...
        ...
        ... % mandatory fields
        ...
        ... % number of simulations for reference data
        'simulationItData', s_numbers.simulationItData, ...
        ... % number of simulations during search
        'simulationItLoop', s_numbers.simulationItLoop, ...
        ... % number of iterations/generated parameter samples
        'iterations', s_numbers.iterations, ...
        ... % size of the epsilon ball
        'epsilon', s_numbers.epsilon, ...
        ... % variance of selection distribution
        'newParameterVariance', s_numbers.newParameterVariance, ...
        ... % *.net file name (without file extension)
        'netFileName', s_strings.netFileName, ...
        ... % measurement time points
        'measurementTimes', s_vectors.measurementTimes, ...
        ... % parameter domains
        'domains', {s_cells.domains}, ...
        ... % measured species
        'measuredVariableNames', {s_cells.measuredVariableNames}, ...
        ... % closure methods used
        'closureMethods', {s_cells.closureMethods}, ...
        ... % true parameters
        'parametersReference', s_objects.parametersReference, ...
        ... % initial parameter estimation
        'parametersInitial', s_objects.parametersInitial, ...
        ...
        ... % fields with optional parameters
        ...
        ... % recompute reference data?
        'computeDataFile', b_computeDataFile, ...
        ... % recompute closure files?
        'computeClosureFiles', b_computeClosureFiles, ...
        ... % recompute simulation script?
        'computeSimulationFile', b_computeSimulationFile, ...
        ... % remove closures after crashing?
        'removeCrashingClosures', b_removeCrashingClosures, ...
        ... % show plots?
        'plot', b_plot, ...
        ... % likelihood computation mode
        'likelihoodMode', n_likelihoodMode, ...
        ... % limit how often new parameters may be rejected
        'rejectionLimit', n_rejectionLimit, ...
        ... % likelihood achieved by closure to be used
        'likelihoodThreshold', n_likelihoodThreshold, ...
        ... % figure placement
        'figurePlacement', n_figurePlacement, ...
        ... % random number generation seed (0 = default, < 0 = none)
        'seed', n_seed, ...
        ... % base name used for several purposes
        'baseName', si_baseName, ...
        ... % simulation function name
        'simulationName', si_simulationName, ...
        ... % figure name
        'figureName', si_figureName, ...
        ... % simulation file
        'simulationFile', si_simulationFile, ...
        ... % logs saving file name
        'diaryName', si_diaryName, ...
        ... % append this suffix to some strings
        'appendName', si_appendName, ...
        ... % data structure for inference of parameter switches
        'dynamicSwitches', v_dynamicSwitches, ...
        ... % measurement times for moments
        'momentTimes', v_momentTimes, ...
        ... % ODE solver function names
        'funNamesOdeSolve', {c_funNamesOdeSolve}, ...
        ... % time-dependent functions
        'timeFunctions', {c_timeFunctions}, ...
        ...
        ... % fields which are set later
        ...
        ... % folder for writing data to
        'dataFolder', [], ...
        ... % net file
        'netFile', [], ...
        ... % reference data file name
        'dataFileName', [], ...
        ... % symbolic parameters
        'parametersSymbolic', [], ...
        ... % simulation indices of measured species
        'simulationIndices', [], ...
        ... % mean indices of measured species
        'meanIndices', [], ...
        ... % variance indices of measured species
        'varianceIndices', [], ...
        ... % covariance index of first and second measured species
        'covarianceIndex', [], ...
        ... % ODE solver functions
        'odeFunctions', [], ...
        ... % ODE solver flags for transposing data
        'odeTransposeData', [], ...
        ... % ODE solver flags for starting at first index
        'odeStartsAtFirstIndex', [], ...
        ... % ODE solver termination bounds
        'odeTerminationBounds', [], ...
        ... % simulation function
        'simulationFunction', [], ...
        ... % moments of the reference data
        'referenceMoments', [], ...
        ... % closure file names
        'closureFileNames', [], ...
        ... % skip simulations during search?
        'skipSimulation', [] ...
    );
end


function [ s_constants ] = setupFiles( s_constants )
%% initializes folder and file names

    si_netFileName = s_constants.netFileName;
    
    % create data folder
    si_dataFolder = ['data', filesep, si_netFileName];
    s_constants.dataFolder = si_dataFolder;
    if (~ exist(si_dataFolder, 'dir'))
        mkdir(si_dataFolder);
    end
    % add it to the path
    c_paths = regexp(path, pathsep, 'Split');
    if (~ ismember(si_dataFolder, c_paths))
        addpath(si_dataFolder);
    end
    
    % net file
    s_constants.netFile = [si_netFileName, '.net'];
    
    % base name
    si_baseName = s_constants.baseName;
    
    % file for reference simulation data
    s_constants.dataFileName = [si_dataFolder, filesep, si_baseName, '_Data'];
    
    % file for simulation function
    si_simulationName = s_constants.simulationName;
    if (isempty(si_simulationName))
        si_simulationName = ['Gillespie', si_baseName];
        s_constants.simulationName = si_simulationName;
    end
    if (isempty(s_constants.simulationFile))
        s_constants.simulationFile = ...
            [si_dataFolder, filesep, si_simulationName, '.m'];
    end
    
    % write diary
    si_diaryPath = s_constants.diaryName;
    if (DEBUG.C_n_storeLogs ~= 0)
        si_diaryPath = [si_dataFolder, filesep, si_diaryPath];
        if (DEBUG.C_n_storeLogs == 2)
            % remove old diary file
            if (exist(si_diaryPath, 'file'))
                delete(si_diaryPath);
            end
        end
        diary(si_diaryPath);
    end
end


function [ s_protocol, s_current, s_constants ] = setupDataStructures( ...
        s_constants, v_netParameters )
%% initialize data structures for the search and the current state
% NOTE: All ParameterVector objects are copied.

    n_params = size(v_netParameters, 1);
    o_paramsInit = s_constants.parametersInitial;
    
    % create a search protocol containing:
    % - all parameters
    % - all closure choices
    % - all events
    s_protocol = repmat(struct(...
            'parameters', [], ...
            'closure', 0, ...
            'isAccepted', 0), ...
        s_constants.iterations, 1);
    s_protocol(1).parameters = ParameterVector(o_paramsInit);
    
    % symbolic parameter names
    v_symbolicParameters = sym(zeros(1, n_params));
    for i = 1 : n_params
        v_symbolicParameters(i) = sym(v_netParameters(i).id);
    end
    s_constants.parametersSymbolic = v_symbolicParameters;
    
    % create a structure containing changing variables:
    % - the current iteration index
    % - the current reference parameters for the stochastic search
    % - the new parameters
    s_current = struct(...
        'index', 1, ...
        'epsilonIndex', 0, ...
        'lastClosureIndex', 0, ...
        'rejectionSeries', 0, ...
        'parametersAccepted', ParameterVector(o_paramsInit), ...
        'parametersCurrent', ParameterVector(o_paramsInit), ...
        'parametersCenter', [], ...
        'parametersBest', ParameterVector(o_paramsInit), ...
        'likelihoodPrevious', -inf, ...
        'likelihoodBest', -inf);
end


function [ s_constants ] = setupDomains( s_constants )
%% initialize data structure for parameter domains

    o_parameters = s_constants.parametersReference;
    n_tMax = s_constants.measurementTimes(end);
    c_domainsInput = s_constants.domains;

    n_totalNo = o_parameters.getTotalParamsNo();
    n_timesNo = o_parameters.getTimesNo();
    n_paramsNo = n_totalNo - n_timesNo;
    if (size(c_domainsInput, 1) < n_totalNo)
        c_domainsNew = cell(n_totalNo, 3);
    else
        assert((size(c_domainsInput, 1) == n_totalNo) && ...
            iscell(c_domainsInput) && ...
            (size(c_domainsInput, 2) == 3), 'Illegal domain input.');
        c_domainsNew = c_domainsInput;
    end
    m_values = [];
    
    for i = 1 : n_totalNo
        for j = 1 : 3
            value = c_domainsNew{i, j};
            if (isempty(value))
                switch j
                    case 1
                        % lower bound
                        value = 0;
                    
                    case 2
                        % upper bound
                        if (i > n_paramsNo)
                            % time parameter
                            value = n_tMax;
                        else
                            % normal parameter
                            value = Inf;
                        end
                    
                    case 3
                        % additional information
                        value = Constants.C_si_domainStandard;
                end
            else
                switch j
                    case {1, 2}
                        if (~ isnumeric(value))
                            error(Utils.wrapError('Illegal domain bound.'));
                        end
                    
                    case 3
                        % domain type
                        switch value
                            case {Constants.C_si_domainStandard, ...
                                    Constants.C_si_domainBoolean}
                                % nothing to do
                            
                            case Constants.C_si_domainConstant
                                % overwrite bounds
                                if (isempty(m_values))
                                    m_values = ...
                                       o_parameters.getAllValuesIncludingTime();
                                end
                                c_domainsNew{i, 1} = m_values(i);
                                c_domainsNew{i, 2} = m_values(i);
                            
                            otherwise
                                error(Utils.wrapError(...
                                    'Illegal domain information.'));
                        end
                end
                continue;
            end
            
            c_domainsNew{i, j} = value;
        end
    end
    
    % overwrite domains
    s_constants.domains = c_domainsNew;
end


function [ s_constants ] = setupSpeciesIndices( s_constants, v_species )
%% finds the indices of mean/variance/covariance for the species to be measured
% 
% data format:
% mu_x1,   ...,    mu_xn, 
% x1x1, x1x2, ..., x1xn,
%       x2x2, ..., x2xn,
%             ...
%                  xnxn
% 
% However, Boolean species have no variances, so they must be skipped.

    n_species = size(v_species, 1);
    c_measuredSpeciesNames = s_constants.measuredVariableNames;
    
    % number of measured species (must be 1 or 2)
    n_measuredSpecies = size(c_measuredSpeciesNames, 2);
    if (n_measuredSpecies == 0)
        error(Utils.wrapError('At least one species must be measured.'));
    end
    if ((s_constants.likelihoodMode == Constants.C_n_likelihoodTwo) && ...
            (n_measuredSpecies ~= 2))
        error(Utils.wrapError(...
            'Only 2 species can be measured in the chosen likelihood mode.'));
    end
    
    % find species indices (= mean indices) for some species index k
    % formula (trivial): k
    v_speciesIndices = zeros(1, n_measuredSpecies);
    for i = 1 : n_measuredSpecies
        si_measuredSpeciesName = c_measuredSpeciesNames{i};
        for j = 1 : n_species
            if (strcmp(si_measuredSpeciesName, v_species(j).id))
                v_speciesIndices(i) = j;
                break;
            end
        end
        assert((v_speciesIndices(i) > 0), ...
            'The specified species could not be found.');
    end
    
    % vector which counts the Boolean variables
    % For Boolean variables there is no variance output from the closures.
    bv_booleans = zeros(1, n_species);
    for i = 1 : n_species
        si_type = v_species(i).type;
        switch (si_type)
            case 'boolean'
                bv_booleans(i) = 1;
            
            case 'stochastic'
                % nothing to do
            
            otherwise
                error(Utils.wrapError(...
                    sprintf('Species of type "%s" are not supported.', ...
                    si_type)));
        end
    end
    
    % find variance indices xkxk for some species index k
    v_varianceIndices = zeros(1, n_measuredSpecies);
    for i = 1 : n_measuredSpecies
        n_idx = v_speciesIndices(i);
        v_varianceIndices(i) = ...
            getVarianceIndex(n_species, n_idx, n_idx, bv_booleans);
    end
    
    % store indices in data structure
    s_constants.measuredSpeciesNumber = n_measuredSpecies;
    s_constants.simulationIndices = v_speciesIndices;
    s_constants.meanIndices = v_speciesIndices;
    s_constants.varianceIndices = v_varianceIndices;
    if (s_constants.likelihoodMode == Constants.C_n_likelihoodTwo)
        % additional covariance for two species measured
        
        % find covariance index xjxk for some species indices j and k.
        n_idx1 = v_speciesIndices(1);
        n_idx2 = v_speciesIndices(2);
        n_covarianceIndex = ...
            getVarianceIndex(n_species, n_idx1, n_idx2, bv_booleans);
        
        s_constants.covarianceIndex = n_covarianceIndex;
    end
    
    % print user information
    fprintf('measuring species %s', c_measuredSpeciesNames{1});
    if (length(c_measuredSpeciesNames) > 1)
        for i = 2 : length(c_measuredSpeciesNames) - 1
            fprintf(', %s', c_measuredSpeciesNames{i});
        end
        if (length(c_measuredSpeciesNames) > 2)
            fprintf(',');
        end
        fprintf(' and %s ', c_measuredSpeciesNames{end});
        switch (s_constants.likelihoodMode)
            case Constants.C_n_likelihoodSingle
                fprintf('assuming pairwise independence\n');

            case Constants.C_n_likelihoodTwo
                fprintf('assuming dependence\n');
        end
    else
        fprintf('\n');
    end
end


function [ n_idx ] = getVarianceIndex( n_species, n_idx1, n_idx2, bv_booleans )
%% computes the (co)variance index given two species indices
% Formula, with n the number of species and j, k the target species indices:
% [n](skip means) + [sum_i=(n-j+1 to n) i](row) + [k - j + 1](column)
% = [n] + [(j - 1) * (2 * n - j + 2) / 2] + [k - j + 1]
% 
% The formula must be modified by subtracting the number of species which are of
% Boolean type. If the input indices are the same (i.e., the output is the
% variance) and the respective species is of Boolean type, the result is a
% special value.

    if (n_idx1 == n_idx2)
        if (bv_booleans(n_idx1))
            error(Utils.wrapError(...
                'Measuring Boolean species is not supported.'));
        else
            % ignore the variance of the respective species (see usage below)
            n_min = n_idx1 - 1;
        end
    else
        n_min = n_idx1;
        if (n_idx2 < n_idx1)
            % make sure the first index is the smaller one
            n_idx1 = n_idx2;
            n_idx2 = n_min;
            n_min = n_idx1;
        end
    end
    
    % standard computation without considering Booleans
    n_idx = (n_species) + ...
        ((n_idx1 - 1) * (2 * n_species - n_idx1 + 2) / 2) + ...
        (n_idx2 - n_idx1 + 1);
    
    % subtract the number of Boolean variances before that index
    for i = 1 : n_min
        if (bv_booleans(i))
            n_idx = n_idx - 1;
        end
    end
end


function [ s_constants ] = setupOdeSolvingFunctions( s_constants )
%% finds the associated function for ODE solving given a keyword

    c_funNameOdeSolve = s_constants.funNamesOdeSolve;
    for i = size(c_funNameOdeSolve, 2) : -1 : 1
        arguments = c_funNameOdeSolve{i};
        if (iscell(arguments))
            assert(size(arguments, 2) == 2, 'There must be two arguments.');
            si_name = arguments{1};
            v_terminationBounds(i) = arguments{2};
        else
            si_name = arguments;
            v_terminationBounds(i) = 0;
        end
        
        switch (si_name)
            case 'matlab'
                fun_func{i} = @momentsMatlab;
                bv_transposeData(i) = false;
                bv_startsAtFirstIndex(i) = false;
            
            case 'sundials'
                fun_func{i} = @momentsSundials;
                bv_transposeData(i) = true;
                bv_startsAtFirstIndex(i) = true;
            
            otherwise
                error('The specified ODE solving name %s is not valid."', ...
                    si_name);
        end
    end
    s_constants.odeFunctions = fun_func;
    s_constants.odeTransposeData = bv_transposeData;
    s_constants.odeStartsAtFirstIndex = bv_startsAtFirstIndex;
    s_constants.odeTerminationBounds = v_terminationBounds;
end


function [ s_constants ] = generateSimulationFile( s_constants, s_net )
%% generates a Gillespie simulation function (file)

    if (s_constants.computeSimulationFile)
        b_generate = true;
    else
        b_generate = (~ exist(s_constants.simulationFile, 'file'));
    end
    
    si_simulation = s_constants.simulationName;
    if (b_generate)
        fprintf('generating simulation file\n');
        generateGillespie(s_constants.simulationFile, si_simulation, s_net);
    else
        fprintf('loading simulation file\n');
    end
    s_constants.simulationFunction = str2func(si_simulation);
end


function [ s_constants ] = setupReferenceMoments( s_constants )
%% simulates for getting reference data or loads it from file; computes moments

    si_dataFileName = s_constants.dataFileName;
    
    if (s_constants.computeDataFile)
        % compute reference data in any case
        b_generate = true;
    else
        % check whether data file exists; in this case skip computation
        b_generate = (~ exist([si_dataFileName, '.mat'], 'file'));
    end
    
    if (b_generate)
        % run reference simulation once and store data in a file
        if (DEBUG.C_b_outputDataFinished)
            fprintf('generating reference data\n');
        end
        [c_referenceData, ~] = getSimulationData(s_constants, ...
            s_constants.parametersReference, s_constants.simulationItData, 0);
        s_referenceDataRaw = struct(...
            'data', {c_referenceData}, ...
            'simulations', s_constants.simulationItData, ...
            'measurementTimes', s_constants.measurementTimes, ...
            'parameterMatrix', ...
                s_constants.parametersReference.getEventValueMatrix() ...
            ); %#ok<NASGU>
        save(si_dataFileName, '-struct', 's_referenceDataRaw');
    else
        % load data from file
        if (DEBUG.C_b_outputDataFinished)
            fprintf('loading reference data\n');
        end
        s_referenceDataRaw = load(si_dataFileName);
        c_referenceData = s_referenceDataRaw.data;
        
        % sanity checks that input settings and data coincide
        
        % number of simulations (test complete)
        if (s_referenceDataRaw.simulations ~= s_constants.simulationItData)
            error(Utils.wrapError(...
                ['The number of simulations specified (', ...
                num2str(s_constants.simulationItData), ...
                ') differs from when the stored data was created (', ...
                num2str(s_referenceDataRaw.simulations), ...
                '). Please run the script again with the correct arguments.']));
        end
        % parameters (test complete up to class type)
        m_parameterMatrix = ...
            s_constants.parametersReference.getEventValueMatrix();
        if (~ isequal(s_referenceDataRaw.parameterMatrix, m_parameterMatrix))
            fprintf('<strong>old parameters:\n');
            disp(s_referenceDataRaw.parameterMatrix);
            fprintf('new parameters:\n');
            disp(m_parameterMatrix);
            fprintf('</strong>');
            error(Utils.wrapError(...
                ['The parameters specified differ from when the stored ', ...
                'data was created. Please run the script again with the ', ...
                'correct arguments.']));
        end
        % measurement times (test complete)
        if (~ isequal(s_referenceDataRaw.measurementTimes, ...
                s_constants.measurementTimes))
            error(Utils.wrapError(...
                ['The measurement times specified\n (', ...
                sprintf('%d ', s_constants.measurementTimes), ...
                ')\ndiffer from when the stored data was created\n (', ...
                sprintf('%d ', s_referenceDataRaw.measurementTimes), ...
               ').\nPlease run the script again with the correct arguments.']));
        end
    end
    [s_referenceMoments, ~] = getSimulationMoments(s_constants, ...
            s_constants.parametersReference, s_constants.simulationItData, ...
            0, c_referenceData);
    
    if (DEBUG.C_b_outputReferenceMeanZero)
        for i = 1 : size(s_referenceMoments.mean, 1)
            for j = 1 : size(s_referenceMoments.mean, 2)
                if (s_referenceMoments.mean(i, j) == 0)
                    Utils.warn('Mean = 0 in the data.\n');
                    break;
                end
            end
        end
    end
    
    s_constants.referenceMoments = s_referenceMoments;
end


function [ s_constants, s_closures ] = setupClosures( s_constants, s_net )
%% computes the closures or load them from previous computation

    % approximations (closures) available
    n_closures = size(s_constants.closureMethods, 1);
    assert((n_closures > 0), 'There must be at least one closure method.');
    c_closureFileNames = cell(n_closures, 1);
    
    s_closures = struct(...
        'number', n_closures, ...
        'name', {cell(n_closures, 1)}, ...
        'degree', zeros(n_closures, 1), ...
        'mdyns', {cell(n_closures, 1)}, ...
        'crashes', []);
    
    % check whether there exists a time parameter
    n_timeFunctions = length(s_constants.timeFunctions);
    n_first = length(s_constants.parametersReference.c_params) + 1;
    c_funcNames = {s_net.parameter(n_first : end).id};
    
    % extra handling for derivative matching closure
    % 0 = not seen any dm closure yet; 1 = already computed; -1 = deactivated
    b_useDmInit = 0;
    
    c_closureMethods = s_constants.closureMethods;
    for i = 1 : n_closures
        % compute the closed moment approximation
        % also create a function file that can be passed to an ODE solver
        % (also compute the open moment system, but we do not need it afterward)
        si_closureName = c_closureMethods{i, 1};
        n_maxDeg = c_closureMethods{i, 2};
        si_closureFileName = ...
            [s_constants.baseName, '_', si_closureName, '_', num2str(n_maxDeg)];
        si_closurePathName = ...
            [s_constants.dataFolder, filesep, si_closureFileName];
        si_closureFile = [si_closurePathName, '.m'];
        si_closureDataPathName = [si_closurePathName, '_Data'];
        c_closureFileNames{i} = si_closureFileName;

        if (s_constants.computeClosureFiles)
            % compute closure in any case
            b_load = false;
        else
            % check whether closure file exists; in this case skip computation
            b_load = (exist(si_closureFile, 'file') && ...
                    exist([si_closureDataPathName, '.mat'], 'file'));
        end

        % compute closure
        if (b_load)
            % load closure from file
            if (DEBUG.C_b_outputClosureLoad)
                fprintf('loading closure %s of degree %d\n', si_closureName, ...
                    n_maxDeg);
            end
            s_mdyn = load(si_closureDataPathName);
        else
            % construct new file
            if ((b_useDmInit ~= -1) && strcmp(si_closureName, 'dm'))
                % repair initial condition for derivative matching closure
                if (b_useDmInit == 0)
                    % compute new initial condition (replace 0 by 0.01)
                    s_netDm = s_net;
                    for j = 1 : length(s_net.species)
                        if (s_net.species(j).initialAmount == 0)
                            b_useDmInit = 1;
                            s_netDm.species(j).initialAmount = 0.01;
                        end
                    end
                    if (b_useDmInit == 0)
                        % no initial condition repair, deactivate
                        b_useDmInit = -1;
                    end
                end
                
                if (b_useDmInit == 1)
                    % swap net structure
                    s_netTmp = s_net;
                    s_net = s_netDm;
                end
            end
            
            if (size(c_closureMethods{i}, 2) == 3)
                % additional volume argument for linear noise closure
                s_mdyn = closureDynamics(s_net, n_maxDeg, ...
                    si_closureName, si_closureFileName, ...
                    s_constants.parametersSymbolic, c_closureMethods{i, 3});
            else
                % no volume argument
                s_mdyn = closureDynamics(s_net, n_maxDeg, ...
                    si_closureName, si_closureFileName, ...
                    s_constants.parametersSymbolic);
            end
            % TODO(T) we could compute the open system only once for each
            % closure degree; sort closure methods by maximum degree then
            
            if ((b_useDmInit == 1) && strcmp(si_closureName, 'dm'))
                % swap net structure again
                s_net = s_netTmp;
            end
            
            % modify file for time parameter
            if (n_timeFunctions > 0)
                si_origFileName = [si_closureFileName, '.m'];
                f_origFile = fopen(si_origFileName, 'r');
                f_outFile = fopen([si_closurePathName, '.m'], 'w');
                
                n_line = 0;
                while (~ feof(f_origFile))
                    si_line = fgetl(f_origFile);
                    n_line = n_line + 1;
                    
                    if (n_line == 1)
                        % insert additional time parameter
                        si_line = strrep(si_line, ')', ',t)\n');
                    elseif (strcmp(si_line, '  % compute function'))
                        % insert additional assignment for time function
                        for j = 1 : n_timeFunctions
                            si_funcName = c_funcNames{j};
                            si_line = sprintf('  %s = %s(t);\n', ...
                                si_funcName, si_funcName);
                            fprintf(f_outFile, si_line, si_line);
                        end
                        continue;
                    else
                        % else replace percent signs and add line delimiter
                        si_line = [strrep(si_line, '%', '%%'), '\n'];
                    end
                    % append line to output file
                    fprintf(f_outFile, si_line);
                end
                
                fclose(f_outFile);
                fclose(f_origFile);
                
                % remove old file
                delete(si_origFileName);
            else
                % move closure file to data folder
                movefile([si_closureFileName, '.m'], si_closureFile);
            end
            
            % generate additional data file
            save(si_closureDataPathName, '-struct', 's_mdyn');
        end
        
        s_closures.name{i} = si_closureName;
        s_closures.degree(i) = n_maxDeg;
        s_closures.mdyns{i} = s_mdyn;
    end
    s_constants.closureFileNames = c_closureFileNames;
    
    % skip simulation in single-closure case
    s_constants.skipSimulation = (s_closures.number <= 1);
    
    % extra moment times for plotting
    v_measurementTimes = s_constants.measurementTimes;
    if (isempty(s_constants.momentTimes))
        % use the same times as for simulation
        s_constants.momentTimes = v_measurementTimes;
        s_constants.momentsFilterIndices = [];
    else
        % use a superset of indices (compared to for simulation)
        % => add a filter at which indices the measurements should take place
        v_momentsFilteredIndices = ...
            zeros(1, length(v_measurementTimes));
        n_currentMomentTimeIndex = 0;
        n_currentTimeIndex = 1;
        n_currentTime = v_measurementTimes(n_currentTimeIndex);
        for n_currentMomentTime = s_constants.momentTimes
            n_currentMomentTimeIndex = n_currentMomentTimeIndex + 1;
            if (n_currentMomentTime == n_currentTime)
                v_momentsFilteredIndices(n_currentTimeIndex) = ...
                    n_currentMomentTimeIndex;
                n_currentTimeIndex = n_currentTimeIndex + 1;
                if (n_currentTimeIndex > length(v_measurementTimes))
                    break;
                else
                    n_currentTime = v_measurementTimes(n_currentTimeIndex);
                end
            elseif (n_currentMomentTime > n_currentTime)
                error(Utils.wrapError(...
                    'Missing time point for moment measurement.'));
            end
        end
        s_constants.momentsFilterIndices = v_momentsFilteredIndices;
    end
end


function [ s_plotSetup ] = setupPlotStructure( )
%% initializes preliminary data structure for plotting

    s_plotSetup = struct(...
        'parameters', true, ...
        'moments', true, ...
        'momentsHidden', true, ...
        'histogram', false ...
        );
end


function [ s_constants ] = setupDynamicSwitches( s_constants )
%% initializes data structure for dynamic switching

    v_dynamicSwitches = s_constants.dynamicSwitches;
    if (isempty(v_dynamicSwitches))
        % no dynamic switching inference
        v_dynamicSwitches = [0, 1, 1, -inf];
        s_constants.dynamicSwitches = v_dynamicSwitches;
    else
        assert((length(v_dynamicSwitches) == 4), ...
            'The vector must contain 4 elements');
        assert(((v_dynamicSwitches(1) > 0) && ...
            (v_dynamicSwitches(1) < ...
                s_constants.parametersReference.getParamsNo()) && ...
            (v_dynamicSwitches(2) <= v_dynamicSwitches(3))), ...
            'The vector for parameter switch inference contains invalid data.');
    end
end


function [ s_statistics ] = setupStatistics( timer_refTotal )
%% initializes statistics structure

    s_statistics = struct(...
        'timeRefTotal', timer_refTotal, ...
        'timeRefSearch', -1, ...
        'timeRefEpsilonBall', -1, ...
        'timeInitialization', toc(timer_refTotal), ...
        'timeSimulation', 0, ...
        'timeOde', 0, ...
        'timeMcmc', 0, ...
        'timePlot', 0, ...
        'timeOdeCrashes', 0);
end