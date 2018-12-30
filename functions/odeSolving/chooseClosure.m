function [ s_protocol, s_current, s_closures, s_statistics ] = ...
    chooseClosure( s_current, s_protocol, s_constants, ...
        s_closures, s_plot, s_statistics )
%% CHOOSECLOSURE Evaluates the closure methods and chooses one to use

    if (DEBUG.C_b_outputReevaluation)
        if (s_current.index > 1)
            n_timeDiff = toc(s_statistics.timeRefEpsilonBall);
            fprintf(...
                'leaving epsilon-ball after %d steps and %.2f seconds\n', ...
                s_current.index - s_current.epsilonIndex, n_timeDiff);
        end
    end
    
    % plot mean and variance for data
    if (s_plot.isPlot)
        n_timeDiff = plotMeanVariance(s_plot, s_constants, ...
            s_constants.referenceMoments, Constants.C_n_plotData);
        s_statistics.timePlot = s_statistics.timePlot + n_timeDiff;
    end
    
    % only run simulation when told so
    if (s_constants.skipSimulation)
        s_sampleMoments = [];
    else
        % run simulation (considered "exact" solution)
        [s_sampleMoments, s_statistics.timeSimulation] = ...
            getSimulationMoments(s_constants, s_current.parametersCurrent, ...
            s_constants.simulationItLoop, s_statistics.timeSimulation, []);
        
        % plot mean and variance for simulation
        if (s_plot.isPlot)
            n_timeDiff = plotMeanVariance(s_plot, s_constants, ...
                s_sampleMoments, Constants.C_n_plotSimulation);
            s_statistics.timePlot = s_statistics.timePlot + n_timeDiff;
        end
    end
    
    % run closure methods (ODE approximations)
    s_closures.likelihoods = zeros(s_closures.number, 1);
    s_closures.crashes = [];
    for j = 1 : s_closures.number
        % compute likelihood wrt. exact solution of current parameters
        [s_closures, s_statistics, n_likelihood, s_odeMoments, ...
            b_crashed] = runClosure(s_current, s_constants, s_closures, ...
            s_statistics, s_sampleMoments, j);
        if (b_crashed)
            continue;
        else
            s_closures.likelihoods(j) = n_likelihood;
        end
        
        % plot mean and variance for closures
        if (s_plot.isPlot)
            n_timeDiff = plotMeanVariance(s_plot, s_constants, s_odeMoments, j);
            s_statistics.timePlot = s_statistics.timePlot + n_timeDiff;
        end
    end
    
    % fix crashes by removing closures
    if (~ isempty(s_closures.crashes))
        [s_constants, s_closures] = fixClosureCrashes(s_constants, s_closures);
        
        assert((length(s_closures.likelihoods) == s_closures.number), ...
            'There should be as many likelihoods as closures remaining.');
        if (s_closures.number == length(s_closures.crashes))
            % trigger termination when no closure worked
            s_current.closureIndex = 0;
            s_closures.number = 0;
            return;
        end
    else
        s_closures.idxOffset = zeros(length(s_closures.name), 1);
    end
    
    % choose the best closure method
    [n_closureIndex, s_current] = findBestLikelihood(s_current, s_closures);
    s_protocol(s_current.index).closure = n_closureIndex;
    
    % plot mean and variance for chosen closure again in bold
    if (s_plot.isPlot)
        n_closureIndexModuloCrashes = ...
            n_closureIndex - s_closures.idxOffset(n_closureIndex);
        n_timeDiff = plotMeanVariance(s_plot, s_constants, [], ...
            -n_closureIndexModuloCrashes, DEBUG.C_b_showLegend);
        s_statistics.timePlot = s_statistics.timePlot + n_timeDiff;
    end
    
    % TODO(T) We might want to ignore parameters where no ODE solution is
    % good and fall back to simulation
%     if (likelihoodClosures(idx) <= likelihoodThreshold)
%         % re-sample if even the best likelihood is not good enough
%         % f.i. if all ODE approxmiations crashed
%         fprintf('ODE likelihood is too small!\n');
%         continue;
%     end
    
    s_current.closureIndex = n_closureIndex;
    s_current.parametersCenter = s_current.parametersCurrent;
    s_current.rejectionSeries = 0;
    s_statistics.timeRefEpsilonBall = tic;
    s_current.epsilonIndex = s_current.index;
end
%% --- end of main function ---


%% --- helper functions ---

function [ n_idx, s_current ] = findBestLikelihood( s_current, s_closures )
%% finds the best likelihood

    v_likelihoodClosures = s_closures.likelihoods;
    n_idx = maxIndex(v_likelihoodClosures);
    
    if (DEBUG.C_b_outputAllClosureLikelihoods)
        fprintf('likelihoods: ');
        for DBG = 1 : length(v_likelihoodClosures)
            if (DBG > 1)
                fprintf(', ');
            end
            if (DBG == n_idx)
                fprintf('<strong>');
            end
            fprintf('%.3f', v_likelihoodClosures(DBG))
            if (DBG == n_idx)
                fprintf('</strong>');
            end
        end
        fprintf('\n');
    end
    
    n_lastClosureIdx = s_current.lastClosureIndex;
    if (n_idx ~= n_lastClosureIdx)
        % new best closure
        s_current.lastClosureIndex = n_idx;
        
        if (DEBUG.C_b_outputClosureChoice)
            if (n_lastClosureIdx == 0)
                fprintf('starting with closure %d (%s deg. %d)\n', ...
                    n_idx, s_closures.name{n_idx}, s_closures.degree(n_idx));
            else
                fprintf(...
    'switching from closure %d (%s deg. %d) to closure %d (%s deg. %d)\n', ...
                    n_lastClosureIdx, s_closures.name{n_lastClosureIdx}, ...
                    s_closures.degree(n_lastClosureIdx), ...
                    n_idx, s_closures.name{n_idx}, s_closures.degree(n_idx));
            end
        end
    else
        % use previous closure again
        if (DEBUG.C_b_outputClosureChoice)
            fprintf('keeping closure %d (%s deg. %d)\n', ...
                n_idx, s_closures.name{n_idx}, ...
                s_closures.degree(n_idx));
        end
    end
end