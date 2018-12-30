function generateGillespie( si_filePath, si_name, s_net )
%% GENERATEGILLESPIE Generates a Gillespie simulation function file
% NOTE: To avoid potential errors if species and/or parameters have a name
% clash with each other and/or variables used by the function itself,
% underscores and prefixes are added automatically.
%
% Input arguments:
% si_path: string which is the path to some folder
% si_name: string which is the name of the function/file
% s_net: structure created by calling readNet() for some net file

    % flag for adding debug output
    b_addDebugOutput = true;
    
    % create and open file
    fid = fopen(si_filePath, 'w');
    
    % common lengths
    n_species = length(s_net.species);
    n_parameters = length(s_net.parameter);
    n_reactions = length(s_net.reaction);
    
    % dictionary for variable names
    s_dict = getDictionary();
    
    % - header -
    fprintf(fid, ...
        'function [ %s ] = %s( %s, %s, %s )\n\n', ...
            s_dict.outData, si_name, s_dict.inParams, s_dict.inIt, ...
            s_dict.inConstants);
    
    % - initialization -
    % time
    fprintf(fid, '%% time constants\n');
    fprintf(fid, '%s = %s.measurementTimes;\n', s_dict.tGrid, ...
        s_dict.inConstants);
    fprintf(fid, '%s = %s(end);\n', s_dict.tMax, s_dict.tGrid);
    fprintf(fid, '%s = (%s(1) == 0);\n', s_dict.tMinZero, s_dict.tGrid);
    fprintf(fid, '%s = length(%s);\n', s_dict.tGridLength, s_dict.tGrid);
    fprintf(fid, '%s = %s.timeFunctions;\n\n', s_dict.timeFunc, ...
        s_dict.inConstants);
    
    % stochiometry matrix
    fprintf(fid, '%% stoichiometry matrix\n');
    fprintf(fid, '%s = ...\n[', s_dict.stochMat);
    for i = 1 : n_reactions
        v_row = getStochiometryRow(s_net, i, n_species);
        for j = 1 : n_species
            n_value = v_row(1, j);
            if (n_value < 0)
                si_tmp = sprintf('\t%d', n_value);
            else
                si_tmp = sprintf('\t %d', n_value);
            end
            if (j < n_species)
                si_tmp = sprintf('%s,', si_tmp);
            elseif (i < n_reactions)
                si_tmp = sprintf('%s;\n', si_tmp);
            end
            fprintf(fid, si_tmp);
        end
    end
    fprintf(fid, '\t]'';\n\n');
    % species
    fprintf(fid, '%% species\n');
    fprintf(fid, '%s = cell(1, %d);\n', s_dict.outData, n_species);
    for i = 1 : n_species
        fprintf(fid, '%s{%d} = zeros(%s, %s);\n', s_dict.outData, i, ...
            s_dict.inIt, s_dict.tGridLength);
    end
    fprintf(fid, '\n');
    % parameters
    fprintf(fid, '%% reaction parameters\n');
    fprintf(fid, '[%s, ~] = %s.getEventValueMatrix();\n', ...
        s_dict.paramsConverted, s_dict.inParams);
    for i = 1 : n_parameters
        fprintf(fid, '%s%s = %s(%d, 1);\n', s_dict.parameterPrefix, ...
            s_net.parameter(i).id, s_dict.paramsConverted, i);
    end
    fprintf(fid, '\n');
    % initial number of species
    fprintf(fid, '%% initial number of species\n');
    fprintf(fid, '%s = [', s_dict.specInit);
    for i = 1 : n_species
        si_tmp = char(s_net.species(i).initialAmount);
        if (i < n_species)
            si_tmp = sprintf('%s; ', si_tmp);
        end
        fprintf(fid, si_tmp);
    end
    fprintf(fid, '];\n\n');
    
    % - for loop -
    fprintf(fid, '%% loop\n');
    fprintf(fid, 'for %s = 1 : %s\n', s_dict.it, s_dict.inIt);
    if (b_addDebugOutput)
        % debug output
        fprintf(fid, '\tif (DEBUG.C_n_outputSimIt > 0)\n');
        fprintf(fid, '\t\tif (mod(%s, DEBUG.C_n_outputSimIt) == 0)\n', ...
            s_dict.it);
        fprintf(fid, '\t\t\tfprintf(''iteration %%d\\n'', %s);\n', s_dict.it);
        fprintf(fid, '\t\tend\n\tend\n\t\n');
    end
    % time reset
    fprintf(fid, '\t%% initialize time variables\n');
    fprintf(fid, '\t%s = %s(1);\n', s_dict.T, s_dict.tGrid);
    fprintf(fid, '\t%s = 0;\n', s_dict.t);
    fprintf(fid, '\t%s = 1;\n\t\n', s_dict.indexT);
    % reset species data structure
    fprintf(fid, '\t%% initialize species\n');
    fprintf(fid, '\t%s = zeros(%d, %s);\n', s_dict.specMatrix, n_species, ...
        s_dict.tGridLength);
    fprintf(fid, '\t%s = %s;\n\t\n', s_dict.specCurrent, s_dict.specInit);
    % begin while loop
    fprintf(fid, '\twhile (true)\n');
    % copy reactions here
    fprintf(fid, '\t\t%% reaction propensities\n');
    for i = 1 : n_reactions
        si_tmp = getReactionSubstitution(s_net, i, n_species, n_parameters, ...
            s_dict);
        fprintf(fid, '\t\t%s%i = %s;\n', s_dict.reactionPrefix, i, si_tmp);
    end
    fprintf(fid, '\t\t%s = ', s_dict.reactionSum);
    for i = 1 : n_reactions
        si_tmp = sprintf('%s%d', s_dict.reactionPrefix, i);
        if (i < n_reactions)
            si_tmp = sprintf('%s + ', si_tmp);
        end
        fprintf(fid, si_tmp);
    end
    fprintf(fid, ';\n\t\t\n');
    % fire time
    fprintf(fid, '\t\t%% simulate the waiting time to the next reaction\n');
    fprintf(fid, '\t\tif (%s == 0)\n', s_dict.reactionSum);
    fprintf(fid, '\t\t\tbreak;\n');
    fprintf(fid, '\t\tend\n');
    fprintf(fid, '\t\t%s = %s + (-log(rand()) / %s);\n', ...
        s_dict.t, s_dict.t, s_dict.reactionSum);
    fprintf(fid, '\t\tif (%s >= %s)\n', s_dict.t, s_dict.tMax);
    fprintf(fid, '\t\t\tbreak;\n');
    fprintf(fid, '\t\tend\n\t\t\n');
    % current state
    fprintf(fid, ['\t\t%% save the trajectory state before the next time ', ...
        'grid point is crossed\n']);
    fprintf(fid, '\t\twhile (%s <= %s)\n', s_dict.T, s_dict.t);
    fprintf(fid, '\t\t\t%s(:, %s) = %s;\n', s_dict.specMatrix, ...
        s_dict.indexT, s_dict.specCurrent);
    fprintf(fid, '\t\t\t%s = %s + 1;\n', s_dict.indexT, s_dict.indexT);
    fprintf(fid, '\t\t\t%s = %s(%s);\n', s_dict.T, s_dict.tGrid, s_dict.indexT);
    fprintf(fid, '\t\tend\n\t\t\n');
    % next simulation
    fprintf(fid, '\t\t%% simulate which reaction fires next\n');
    fprintf(fid, '\t\t%s = [', s_dict.p);
    for i = 1 : n_reactions
        si_tmp = sprintf('%s%d / %s', s_dict.reactionPrefix, i, ...
            s_dict.reactionSum);
        if (i < n_reactions)
            si_tmp = sprintf('%s, ...\n\t\t\t', si_tmp);
        end
        fprintf(fid, si_tmp);
    end
    fprintf(fid, '];\n');
    fprintf(fid, '\t\t%s = [0 cumsum(%s)];\n', s_dict.cumprob, s_dict.p);
    fprintf(fid, '\t\t%s = find(rand() > %s, 1, ''last'');\n\t\t\n', ...
        s_dict.indexReaction, s_dict.cumprob);
    fprintf(fid, '\t\t');
    % execute reaction
    fprintf(fid, '%% update the state using the stoichiometry matrix\n');
    fprintf(fid, '\t\t%s = %s + %s(:, %s);\n', s_dict.specCurrent, ...
        s_dict.specCurrent, s_dict.stochMat, s_dict.indexReaction);
    % end while loop
    fprintf(fid, '\tend\n\t\n');
    % finish writing data for last indices
    fprintf(fid, '\t%% write remaining data\n');
    fprintf(fid, '\twhile (%s < %s)\n', s_dict.T, s_dict.tMax);
    fprintf(fid, '\t\t%s(:, %s) = %s;\n', s_dict.specMatrix, ...
        s_dict.indexT, s_dict.specCurrent);
    fprintf(fid, '\t\t%s = %s + 1;\n', s_dict.indexT, s_dict.indexT);
    fprintf(fid, '\t\t%s = %s(%s);\n', s_dict.T, s_dict.tGrid, s_dict.indexT);
    fprintf(fid, '\tend\n');
    fprintf(fid, '\t%s(:, %s) = %s;\n\t\n', s_dict.specMatrix, ...
        s_dict.tGridLength, s_dict.specCurrent);
    % store species
    fprintf(fid, '\t%% overwrite initial values\n');
    fprintf(fid, '\tif (%s)\n', s_dict.tMinZero);
    fprintf(fid, '\t\t%s(:, 1) = %s;\n', s_dict.specMatrix, s_dict.specInit);
    fprintf(fid, '\tend\n\t\n');
    fprintf(fid, '\t%% save the trajectory\n');
    for i = 1 : n_species
        fprintf(fid, '\t%s{1, %d}(%s, :) = %s(%d, :);\n', s_dict.outData, i, ...
            s_dict.it, s_dict.specMatrix, i);
    end
    fprintf(fid, 'end\n');
    
    % - finished file writing -
    fclose(fid);
end

function [ s_dict ] = getDictionary( )
%% creates a dictionary with all variable names used

    s_dict = struct(...
        'outData', 'c_data', ...
        'inParams', 'o_params', ...
        'inIt', 'n_iterations', ...
        'inConstants', 's_constants', ...
        'timeFunc', 'c_timeFunctions', ...
        'it', 'i_', ...
        't', 't_', ...
        'T', 'T_', ...
        'indexT', 'indexT_', ...
        'p', 'x_', ...
        'cumprob', 'x_', ...
        'indexReaction', 'x_', ...
        'tMax', 'tMax_', ...
        'tGrid', 'v_timegrid', ...
        'tGridLength', 'tGridLength_', ...
        'paramsConverted', 'parameters_', ...
        'stochMat', 'stoch_', ...
        'specInit', 'spec_init', ...
        'specMatrix', 'spec_matrix', ...
        'specCurrent', 'spec_cur', ...
        'reactionSum', 'sum_', ...
        'tMinZero', 'b_minEqualsZero', ...
        'speciesPrefix', 's_', ...
        'parameterPrefix', 'p_', ...
        'reactionPrefix', 'h_' ...
        );
end

function [ v_row ] = getStochiometryRow( s_net, n_row, n_colMax )
%% computes the i-th row of the stochiometry matrix from the net data structure

    s_reaction = s_net.reaction(n_row);
    c_species = s_net.species;
    v_row = zeros(1, n_colMax);
    for i = 1 : length(s_reaction.old)
        % find species index
        si_old = s_reaction.old{i};
        n_speciesIdx = 0;
        for j = 1 : length(c_species)
            if (strcmp(si_old, c_species(j).id))
                n_speciesIdx = j;
                break;
            end
        end
        assert(n_speciesIdx > 0, 'The species must exist.');
        
        % extract suffix
        si_new = s_reaction.new{i};
        si_suffix = si_new((length(si_old) + 1) : end);
        
        % remove space characters
        si_suffix = regexprep(si_suffix, '\s', '');
        
        % check whether it is increasing or decreasing
        assert(~ isempty(si_suffix), ...
            'The reaction must specify a change to the species.');
        
        % addition or subtraction
        if (strcmp(si_suffix(1), '+'))
            n_factor = 1;
        elseif (strcmp(si_suffix(1), '-'))
            n_factor = -1;
        else
            assert(false, 'The reaction must contain a ''+'' or ''-'' sign.');
        end
        
        % amount of change
        v_row(1, n_speciesIdx) = n_factor * (str2num(si_suffix(2 : end))); ...
            %#ok<ST2NM>
    end
end

function [ si_out ] = getReactionSubstitution( ...
        s_net, n_idx, n_species, n_parameters, s_dict )
%% substitutes the species and parameter names in the reaction string
% NOTE: To simulate parallel substitution, we have a two-level substitution

    % replace species and parameters by temporary placeholders
    si_out = s_net.reaction(n_idx).intensity;
    for i = 1 : n_species
        si_speciesName = s_net.species(i).id;
        si_out = strrep(si_out, si_speciesName, ...
            sprintf('__&(%d)', i));
    end
    for i = 1 : n_parameters
        si_parameterName = s_net.parameter(i).id;
        si_out = strrep(si_out, si_parameterName, ...
            sprintf('%s%s', s_dict.parameterPrefix, si_parameterName));
    end
    
    % replace temporary placeholders by final strings
    si_out = strrep(si_out, '__&', s_dict.specCurrent);
end