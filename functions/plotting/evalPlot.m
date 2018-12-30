function evalPlot( s_current, s_constants, s_closures, s_plot )
%% EVALPLOT Plots the results

    fprintf('\nplotting results\n');
    
    s_closures.likelihoods = zeros(s_closures.number, 1);
    s_closures.crashes = [];
    
    for i = 1 : 2
        if (i == 1)
            % plot for best parameters found
            o_parameters = s_current.parametersBest;
            
            % rename figure
            set(s_plot.figMeanVar.fig, ...
                'Name', 'mean / variance best');
        else
            % plot for reference parameters
            o_parameters = s_constants.parametersReference;
            
            % use new figures for comparison at reference parameters
            % For this we set the pointer to the other figure.
            set(s_plot.figMeanVarRef.fig, 'Visible', 'on');
            s_plot.figMeanVar = s_plot.figMeanVarRef;
        end
        s_current.parametersCurrent = o_parameters;
        
        % plot mean and variance for data
        plotMeanVariance(s_plot, s_constants, s_constants.referenceMoments, ...
            Constants.C_n_plotData);
        
        % plot mean and variance of simulation for best parameters
        if (i == 1)
            s_sampleMoments = getSimulationMoments(s_constants, ...
                o_parameters, s_constants.simulationItLoop, 0, []);
            plotMeanVariance(s_plot, s_constants, s_sampleMoments, ...
                Constants.C_n_plotSimulation);
        else
            s_sampleMoments = s_constants.referenceMoments;
        end
        
        % no closure computation if there is none
        if (s_closures.number == 0)
            continue;
        end
        
        % run closure methods again
        c_odeMoments = cell(s_closures.number, 1);
        for j = 1 : s_closures.number
            % compute likelihood wrt. exact solution of current parameters
            [s_closures, ~, n_likelihood, s_odeMoments, b_crashed] = ...
                runClosure(s_current, s_constants, s_closures, [], ...
                s_sampleMoments, j);
            if (b_crashed)
                continue;
            else
                s_closures.likelihoods(j) = n_likelihood;
                c_odeMoments{j} = s_odeMoments;
            end
        end
        
        % fix crashes by removing closures
        if (~ isempty(s_closures.crashes))
            [s_constants, s_closures] = fixClosureCrashes(...
                s_constants, s_closures);
            assert((length(s_closures.likelihoods) == s_closures.number), ...
                'There should be as many likelihoods as closures remaining.');
            if ((s_closures.number == 0) || ...
                    ((~ s_constants.removeCrashingClosures) && ...
                    (s_closures.number == length(s_closures.crashes))))
                if (i == 1)
                    attribute = 'reference';
                else
                    attribute = 'best';
                end
                fprintf(...
                    'no working closure exists for the %s parameters...\n', ...
                    attribute);
                continue;
            end
            s_closures.crashes = [];
        end
        
        % choose the best closure method
        n_closureIndex = maxIndex(s_closures.likelihoods);
        
        % plot mean and variance for chosen closure again in bold
        v_closureMoments = c_odeMoments{n_closureIndex};
        if (~ isempty(v_closureMoments))
            plotMeanVariance(s_plot, s_constants, ...
                c_odeMoments{n_closureIndex}, Constants.C_n_plotBestClosure, ...
                DEBUG.C_b_showLegend);
        end
    end
end