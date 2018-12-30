function [ n_time ] = plotMeanVariance( s_plot, s_constants, s_moments, ...
        n_mode, b_showLegend )
%% PLOTMEANVARIANCE Plot mean and variance for data, simulation, and closures
%
% mode:
% 0 <= mode < 1: add plot (special cases: data, simulation, best closure)
% 1 <= mode: add closure plot
% mode < 0: update plot
%
% NOTE: The 'CurrentFigure' trick prevents the figure stealing the window focus.

    timer_all = tic;
    
    set(0, 'CurrentFigure', s_plot.figMeanVar.fig);
    
    if (n_mode >= 0)
        % add a graph to the plot
        
        % determine plot properties
        switch (n_mode)
            case Constants.C_n_plotData % data
                si_color = ['r', DEBUG.C_si_dataPlotMode];
                n_lineWidth = 3;
                si_legendText = 'reference';
                v_timePoints = s_constants.measurementTimes;
            case Constants.C_n_plotSimulation % simulation
                si_color = 'b';
                n_lineWidth = 3;
                si_legendText = 'simulation';
                v_timePoints = s_constants.measurementTimes;
            case Constants.C_n_plotBestClosure % best closure
                si_color = 'k';
                n_lineWidth = 3;
                si_legendText = 'best closure';
                v_timePoints = s_constants.momentTimes;
            otherwise % closures
                si_colors = 'rbgcm';
                si_color = si_colors(1 + mod(n_mode, length(si_colors)));
                n_lineWidth = 1;
                si_legendText = ['closure ', num2str(n_mode)];
                v_timePoints = s_constants.momentTimes;
        end
        
        % plot data
        for i = length(s_constants.simulationIndices) : -1 : 1
            n_idx = 2 * i;
            c_data(:, n_idx) = s_moments.variance(:, i);
            c_data(:, n_idx - 1) = s_moments.mean(:, i);
        end
        if (isfield(s_moments, 'covariance'))
            c_data(:, end + 1) = s_moments.covariance;
        end
        assert(size(c_data, 2) == length(s_plot.figMeanVar.subplots), ...
            'The number of data and sub-plots must be the same.');
        for i = 1 : size(c_data, 2)
            subplotHandle = s_plot.figMeanVar.subplots(i);
            set(s_plot.figMeanVar.fig, 'CurrentAxes', subplotHandle);
            if (n_mode == 0)
                % "hold off" removes title, this command keeps it
                set(subplotHandle, 'NextPlot', 'replacechildren');
            else
                hold on
            end
            plot(v_timePoints, c_data(:, i), si_color, ...
                'LineWidth', n_lineWidth, ...
                'DisplayName', si_legendText);
        end
    else
        % update the graph of the given index
        v_allChildren = get(s_plot.figMeanVar.subplots(1), 'Children');
        n_idx = length(v_allChildren) + n_mode;
        if (~ s_constants.skipSimulation)
            % additional plot of simulation data
            n_idx = n_idx - 1;
        end
        for i = 1 : length(s_plot.figMeanVar.subplots)
            subplotHandle = s_plot.figMeanVar.subplots(i);
            set(s_plot.figMeanVar.fig, 'CurrentAxes', subplotHandle);
            
            % draw closure in black and bold
            v_allChildren = get(subplotHandle, 'Children');
            set(v_allChildren(n_idx), 'Color', 'k', 'LineWidth', 4, ...
                'DisplayName', [get(v_allChildren(n_idx), ...
                    'DisplayName'), '*']);
            
            % draw data and simulation in bold
            set(v_allChildren(end), 'LineWidth', 4);
            if (~ s_constants.skipSimulation)
                set(v_allChildren(end - 1), 'LineWidth', 4);
            end
        end
    end
    
    % draw/update legend if enabled (optional argument)
    if (exist('b_showLegend', 'var') && b_showLegend)
        for i = 1 : length(s_plot.figMeanVar.subplots)
            subplotHandle = s_plot.figMeanVar.subplots(i);
            set(s_plot.figMeanVar.fig, 'CurrentAxes', subplotHandle);
            legend('-DynamicLegend', 'Location', 'best');
        end
    end
    
    % update plots
    drawnow;
    
    n_time = toc(timer_all);
end