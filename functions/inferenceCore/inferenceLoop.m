function [ s_current, s_protocol, s_statistics ] = ...
    inferenceLoop( s_current, s_protocol, s_constants, s_plot, ...
        s_closures, s_statistics )
%INFERENCELOOP Iterate the inference loop of the basic algorithm

% the iteration number is global so it can be seen from other contexts
global G_n_iteration;

% time measuring
timer_refLoop = tic;

for G_n_iteration = 1 : s_constants.iterations
    if (DEBUG.C_b_haltFigure)
        debugInterrupt(2);
    end
    if ((DEBUG.C_n_outputTotalIt > 0) && ...
            (mod(G_n_iteration, DEBUG.C_n_outputTotalIt) == 0))
        fprintf('iteration %d, elapsed time = %.0f seconds\n', ...
            G_n_iteration, toc(s_statistics.timeRefSearch));
    end
    
    s_current.index = G_n_iteration;
    
    if (~ isempty(s_closures.crashes))
        % current closure has crashed, need to reevaluate
        s_closures.crashes = [];
        b_reevaluate = true;
    else
        b_reevaluate = ...
            checkRecomputeClosure(s_current, s_constants, s_closures);
    end
    
    if (b_reevaluate)
        % reevaluate the closure method
        [s_protocol, s_current, s_closures, s_statistics] = ...
            chooseClosure(s_current, s_protocol, s_constants, s_closures, ...
            s_plot, s_statistics);
        
        if (s_closures.number == 0)
            fprintf('no working closures remain, terminating in iteration %d...\n', ...
                G_n_iteration);
            break;
        end
    else
        % new parameters sample (stochastic search)
        [s_current, s_statistics.timeMcmc] = ...
            newParameters(s_current, s_constants, s_statistics.timeMcmc);
    end
    
    % compute likelihood wrt. data
    s_closures.crashes = [];
    [s_closures, s_statistics, n_likelihoodNew, s_odeMoments, ~] = ...
        runClosure(s_current, s_constants, s_closures, s_statistics, ...
        s_constants.referenceMoments, s_current.closureIndex);
    
    % fix crash by removing closure
    if (~ isempty(s_closures.crashes))
        [s_constants, s_closures] = fixClosureCrashes(s_constants, s_closures);
        if (s_closures.number <= 1)
            fprintf(...
                'no working closures remain, terminating in iteration %d...\n', ...
                G_n_iteration);
            break;
        end
        s_protocol(G_n_iteration).parameters = s_current.parametersAccepted;
        continue;
    end
    
    % accept or reject new parameter values
    if (b_reevaluate || isAccept(n_likelihoodNew, s_current))
        % accept new parameter values in two cases:
        % 1) The closure was chosen. Then there was no parameter change.
        % 2) The acceptance probability was big enough.
        
        if (DEBUG.C_b_outputAccepted)
            fprintf('Accepted!\n');
        end
        
        s_current.parametersAccepted = s_current.parametersCurrent;
        s_protocol(G_n_iteration).isAccepted = true;
        s_current.likelihoodPrevious = n_likelihoodNew;
        if (s_current.likelihoodBest < n_likelihoodNew)
            s_current.likelihoodBest = n_likelihoodNew;
            s_current.parametersBest = s_current.parametersCurrent;
        end
        s_current.rejectionSeries = 0;
    else
        % reject new parameter values
        
        if (DEBUG.C_b_outputAccepted)
            fprintf('Rejected!\n');
        end
        
        s_protocol(G_n_iteration).isAccepted = false;
        s_current.rejectionSeries = s_current.rejectionSeries + 1;
    end
    s_protocol(G_n_iteration).parameters = s_current.parametersAccepted;
    
    % plotting
    if (s_plot.isPlot)
        % plot parameters
        if (s_plot.lastPlotIndex + DEBUG.C_n_plotParams <= G_n_iteration)
            timeDiff = plotParameters(s_plot, s_protocol, s_current, ...
                s_constants, false);
            s_plot.lastPlotIndex = G_n_iteration;
            s_statistics.timePlot = s_statistics.timePlot + timeDiff;
        end
        
        % plot mean and variance for data and closure
        if ((DEBUG.C_n_plotMeanVar > 0) && ...
                (mod(G_n_iteration, DEBUG.C_n_plotMeanVar) == 0))
            assert(~ isempty(s_odeMoments), ...
                'The closure must have succeeded here.');
            n_timeDiff = plotMeanVariance(s_plot, s_constants, ...
                s_constants.referenceMoments, Constants.C_n_plotData);
            n_timeDiff = n_timeDiff + plotMeanVariance(s_plot, ...
                s_constants, s_odeMoments, ...
                Constants.C_n_plotBestClosure, DEBUG.C_b_showLegend);
            s_statistics.timePlot = s_statistics.timePlot + n_timeDiff;
        end
    end
end

if (DEBUG.C_b_outputTotalTime)
    fprintf('finished inference loop after %.2f seconds\n', ...
        toc(timer_refLoop));
end

% show remaining plots
if (s_plot.isPlot)
    % if the run has crashed at the last index, ignore it
    b_ignoreLastIndex = (~ isempty(s_closures.crashes));
    
    n_newDelta = G_n_iteration - s_plot.lastPlotIndex - b_ignoreLastIndex;
    
    if (n_newDelta > 0)
        timeDiff = plotParameters(s_plot, s_protocol, s_current, ...
            s_constants, false, n_newDelta, b_ignoreLastIndex);
        s_statistics.timePlot = s_statistics.timePlot + timeDiff;
    end
end

end