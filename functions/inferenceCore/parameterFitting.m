function parameterFitting( ...
        s_flags, s_numbers, s_strings, s_vectors, s_cells, s_objects)
%% PARAMETERFITTING Parameter fitting algorithm
% There are many different versions in which this function can be used.

% ----------------------
% --- initialization ---
% ----------------------

if (DEBUG.C_b_outputSkeleton)
    fprintf('------------------\n- initialization -\n------------------\n');
end

% initialization of common data
[s_constants, s_protocol, s_current, s_statistics, s_closures, ...
    s_plotSetup] = initialize(s_flags, s_numbers, s_strings, s_vectors, ...
    s_cells, s_objects);
% clear input variables to save memory and avoid confusions
clear s_flags s_numbers s_strings s_vectors s_cells s_objects;

% statistics timing
s_statistics.timeRefSearch = tic;

% dynamic switches
v_dynamicSwitches = s_constants.dynamicSwitches;

% other variables
s_currentBest = [];
s_protocolBest = [];

if (DEBUG.C_b_outputSkeleton)
    fprintf('finished initialization after %.2f seconds\n', ...
        s_statistics.timeInitialization);
end

% --------------------------------
% --- parameter inference loop ---
% --------------------------------

% for each specified number of switches apply the inference algorithm
for n_iterationSwitch = v_dynamicSwitches(2) : v_dynamicSwitches(3)
    if (DEBUG.C_b_outputSkeleton)
        fprintf('-----------------\n- starting loop -\n-----------------\n');
    end
    
    if (DEBUG.C_b_outputUpdateSwitch)
        fprintf('turning parameter %d into parameter with %d switches\n', ...
            v_dynamicSwitches(1), n_iterationSwitch);
    end
    
    % update designated switching parameter
    [s_currentNew, s_protocolNew, v_refMap] = updateSwitch(s_current, ...
        s_protocol, v_dynamicSwitches(1), n_iterationSwitch, s_constants.seed);
    
    % new plot struct
    s_plot = setupPlots(s_currentNew, v_refMap, s_constants, s_plotSetup);
    
    % main inference algorithm for current fixed parameters
    [s_currentNew, s_protocolNew, s_statistics] = inferenceLoop(...
        s_currentNew, s_protocolNew, s_constants, s_plot, s_closures, ...
        s_statistics);
    
    % evaluation of current inference results
    [s_currentBest, s_protocolBest, b_continueSwitchInference] = ...
        evaluateInference(s_currentNew, s_currentBest, s_protocolNew, ...
        s_protocolBest, v_dynamicSwitches(4));
    
    if (~ b_continueSwitchInference)
        if (DEBUG.C_b_outputSkeleton)
            fprintf(' Not enough improvement ... terminating\n');
        end
        
        break;
    end
end

% -------------------------------
% --- presentation of results ---
% -------------------------------

% print statistics
if (DEBUG.C_b_outputResults)
    printResults(s_currentBest, s_constants, s_statistics);
end

% plot mean and variance for best results
if (s_plot.isPlot)
    evalPlot(s_currentBest, s_constants, s_closures, s_plot);
end

% save data for the future in a new folder
if (DEBUG.C_b_storeResults)
    s_protocol = s_protocolBest; %#ok<NASGU>
    s_current = s_currentBest; %#ok<NASGU>
    c_fileNames = {'s_protocol', 's_current', 's_constants', 's_plot'};
    si_folderName = ...
        ['results', filesep, s_constants.baseName, s_constants.appendName, '-'];
    n_folderIdx = 1;
    while (exist([si_folderName, num2str(n_folderIdx)], 'dir'))
        n_folderIdx = n_folderIdx + 1;
    end
    si_folderName = [si_folderName, num2str(n_folderIdx)];
    fprintf('storing results in folder %s\n', si_folderName);
    mkdir(si_folderName);
    for i = 1 : length(c_fileNames)
        save([si_folderName, filesep, c_fileNames{i}], c_fileNames{i});
    end

    % save figures
    for i = 1 : 3
        figure(i);
        si_name = sprintf('%s%s%s-%d', si_folderName, filesep, ...
            s_constants.figureName, i);
        storeFigure(si_name);
    end
end

if (DEBUG.C_n_storeLogs ~= 0)
    diary off;
    
    if (DEBUG.C_b_storeResults)
        % move diary file to results folder
        si_diaryName = s_constants.diaryName;
        si_diaryPath = [s_constants.dataFolder, filesep, si_diaryName];
        si_target = [si_folderName, filesep, si_diaryName];
        movefile(si_diaryPath, si_target);
    end
end

end
%% --- end of main function ---


%% --- helper functions ---


function printResults( s_current, s_constants, s_statistics )
%% print statistics

    fprintf('----------\n- result -\n----------\n');
    v_paramsBest = s_current.parametersBest.getAllValuesIncludingTime();
    v_paramsRef = s_constants.parametersReference.getAllValuesIncludingTime();
    minLength = min(length(v_paramsRef), length(v_paramsBest));
    n_totalDifferenceRelative = 0;
    for i = 1 : minLength
        n_best = v_paramsBest(i);
        n_reference = v_paramsRef(i);
        n_difference = abs(n_best - n_reference);
        n_differenceRelative = n_difference / n_reference * 100;
        fprintf('p%d  = %.10f\n', i, n_best);
        fprintf('p%d* = %.10f\n', i, n_reference);
        fprintf('diff: %.2f %%\n', n_differenceRelative);
        n_totalDifferenceRelative = ...
            n_totalDifferenceRelative + n_differenceRelative;
    end
    fprintf('average diff: %.2f %%\n', n_totalDifferenceRelative / minLength);
    % parameters which only occur in one vector
    if (length(v_paramsRef) < length(v_paramsBest))
        v_paramsBigger = v_paramsBest;
    else
        v_paramsBigger = v_paramsRef;
    end
    for i = (minLength + 1) : length(v_paramsBigger)
        fprintf('p%d  = %.10f\n', i, v_paramsBigger(i));
    end
    
    n_timeTotal = toc(s_statistics.timeRefTotal);
    n_timeSearch = toc(s_statistics.timeRefSearch);
    fprintf('\ntotal time: %.2f seconds (~ %.0f minutes), thereof:\n', ...
        n_timeTotal, n_timeTotal / 60);
    fprintf('1) initialization time: %.2f seconds\n', ...
        s_statistics.timeInitialization);
    fprintf('2) search time: %.2f seconds, thereof:\n', n_timeSearch);
    fprintf('2.1) simulation time: %.2f seconds\n', ...
        s_statistics.timeSimulation);
    fprintf('2.2) ODE solving time: %.2f seconds\n', s_statistics.timeOde);
    fprintf('2.3) MCMC time: %.2f seconds\n', s_statistics.timeMcmc);
    fprintf('2.4) ODE solver crashes time: %.2f seconds\n', ...
        s_statistics.timeOdeCrashes);
    fprintf('2.5) plotting time: %.2f seconds\n', s_statistics.timePlot);
    n_timeLost = n_timeSearch - s_statistics.timeSimulation - ...
        s_statistics.timeOde - s_statistics.timeMcmc - ...
        s_statistics.timePlot - s_statistics.timeOdeCrashes;
    fprintf('2.6) untracked time in search: %.2f seconds\n', n_timeLost);
    if (n_timeLost > (n_timeTotal / 10))
        fprintf('<strong>NOTE: This is %.2f percent!</strong>\n', ...
            100 / n_timeTotal * n_timeLost);
    end
end