function mainWrapper( ...
        s_flags, s_numbers, s_strings, s_vectors, s_cells, s_objects )
%MAINWRAPPER Wraps the parameter identification algorithm
% This function enables fast pre-/postprocessing of whole benchmark suites.
% The behavior is controlled by global options set by the caller.
% 
% NOTE: Wrapper code must be written in a form that it can call itself
% recursively. For such recursive calls respective functions must deactivate
% their own trigger option to avoid nontermination.

    % --- global options ---
    
    global G_b_speciesMeasuring;
    global G_b_singleClosureEvaluation;
    global G_b_epsilonEvaluation;
    
    % --- direct options ---
    
    % none so far
    
    % --- recursive options ---
    
    if (G_b_speciesMeasuring)
        wrapSpeciesMeasuring(s_flags, s_numbers, s_strings, s_vectors, ...
            s_cells, s_objects);
        return;
    end
    
    if (G_b_singleClosureEvaluation)
        wrapClosureEvaluation(s_flags, s_numbers, s_strings, ...
            s_vectors, s_cells, s_objects);
        return;
    end
    
    if (G_b_epsilonEvaluation)
        wrapEpsilonEvaluation(s_flags, s_numbers, s_strings, ...
            s_vectors, s_cells, s_objects);
        return;
    end
    
    % --- run parameter identification procedure ---
    % Exception handling is necessary, or else global options will not be reset.
    try
        parameterFitting(s_flags, s_numbers, s_strings, s_vectors, s_cells, ...
            s_objects);
    catch exception
        % print exception
        fprintf(getReport(exception));
        
        % deactivate log writing
        if (DEBUG.C_n_storeLogs ~= 0)
            diary off;
        end
    end
end


function wrapSpeciesMeasuring( ...
        s_flags, s_numbers, s_strings, s_vectors, s_cells, s_objects )
% recursive option: measure species in different ways

    fprintf('<strong>- Benchmark mode: Species measuring -</strong>\n');
    
    % deactivate option for recursive calls
    global G_b_speciesMeasuring;
    G_b_speciesMeasuring = false;
    
    s_net = readNet([s_strings.netFileName, '.net']);
    v_species = s_net.species;
    c_speciesNames = cell(1, length(v_species));
    for i = 1 : length(v_species)
        c_speciesNames{i} = v_species(i).id;
    end
    n_species = length(c_speciesNames);
    bv_booleans = zeros(1, n_species);
    for i = 1 : n_species
        if (strcmp(v_species(i).type, 'boolean'))
            bv_booleans(i) = 1;
        end
    end
    
    % 1) only one species
    s_numbers.likelihoodMode = Constants.C_n_likelihoodSingle;
    for i = 1 : n_species
        if (bv_booleans(i))
            continue;
        end
        s_cells.measuredVariableNames = c_speciesNames(i);
        
        % recursive call
        mainWrapper(s_flags, s_numbers, s_strings, s_vectors, s_cells, ...
            s_objects);
    end
    
    % 2) two species (any combination)
    for i = 1 : n_species
        if (bv_booleans(i))
            continue;
        end
        for j = (i + 1) : n_species
            if (bv_booleans(j))
                continue;
            end
            s_cells.measuredVariableNames = ...
                {c_speciesNames{i}, c_speciesNames{j}};
            
            for k = 1 : 2
                if (k == 1)
                    % 2.1) independently
                    s_numbers.likelihoodMode = Constants.C_n_likelihoodSingle;
                else
                    % 2.2) dependently
                    s_numbers.likelihoodMode = Constants.C_n_likelihoodTwo;
                end
                
                % recursive call
                mainWrapper(s_flags, s_numbers, s_strings, s_vectors, ...
                    s_cells, s_objects);
            end
        end
    end
    
    % 3) all species (if more than two)
    if (n_species > 2)
        s_cells.measuredVariableNames = c_speciesNames;
        for i = n_species : -1 : 1
            if (bv_booleans(i))
                s_cells.measuredVariableNames(end) = [];
            end
        end
        if (length(s_cells.measuredVariableNames) > 2)
            s_numbers.likelihoodMode = Constants.C_n_likelihoodSingle;
            % recursive call
            mainWrapper(s_flags, s_numbers, s_strings, s_vectors, s_cells, ...
                s_objects);
        end
    end
    
    % activate option again
    G_b_speciesMeasuring = true;
end


function wrapClosureEvaluation( ...
        s_flags, s_numbers, s_strings, s_vectors, s_cells, s_objects )
%% recursive option: evaluates each closure on its own

    fprintf('<strong>- Benchmark mode: Closure evaluation -</strong>\n');
    
    % deactivate option for recursive calls
    global G_b_singleClosureEvaluation;
    G_b_singleClosureEvaluation = false; %#ok<NASGU>
    
    if (isfield(s_strings, 'appendName'))
        si_appendName = s_strings.appendName;
    else
        si_appendName = '';
    end
    c_closureMethods = s_cells.closureMethods;
    n_closureMethods = size(c_closureMethods, 1);
    for i = 1 : n_closureMethods
        % use a single closure
        c_currentClosure = c_closureMethods(i, :);
        s_cells.closureMethods = c_currentClosure;
        s_strings.appendName = [si_appendName, '-', c_currentClosure{1}, ...
            num2str(c_currentClosure{2})];
        
        % recursive call
        mainWrapper(...
            s_flags, s_numbers, s_strings, s_vectors, s_cells, s_objects);
    end
    
    % use all closures now
    s_cells.closureMethods = c_closureMethods;
    s_strings.appendName = [si_appendName, '-adaptive'];
    
    % recursive call
    mainWrapper(...
        s_flags, s_numbers, s_strings, s_vectors, s_cells, s_objects);
    
    % activate option again
    G_b_singleClosureEvaluation = true;
end


function wrapEpsilonEvaluation( ...
        s_flags, s_numbers, s_strings, s_vectors, s_cells, s_objects )
%% recursive option: evaluates several epsilon values

    fprintf('<strong>- Benchmark mode: epsilon evaluation -</strong>\n');
    
    % deactivate option for recursive calls
    global G_b_epsilonEvaluation;
    global G_n_epsilonEvaluation;
    G_b_epsilonEvaluation = false;
    
    if (isfield(s_strings, 'appendName'))
        si_appendName = s_strings.appendName;
    else
        si_appendName = '';
    end
    
    if (isempty(G_n_epsilonEvaluation))
        error(Utils.wrapError('No epsilon values specified.'));
        return;
    end
    
    epsilonTmp = G_n_epsilonEvaluation;
    for i = 1 : length(epsilonTmp)
        % extract next epsilon value
        n_epsilon = epsilonTmp(i);
        
        s_numbers.epsilon = n_epsilon;
        
        s_strings.appendName = [si_appendName, '-epsilon', num2str(n_epsilon)];
        fprintf('epsilon = %d\n', n_epsilon);
        
        % recursive call
        mainWrapper(...
            s_flags, s_numbers, s_strings, s_vectors, s_cells, s_objects);
    end
    
    % activate option again
    G_b_epsilonEvaluation = true;
end