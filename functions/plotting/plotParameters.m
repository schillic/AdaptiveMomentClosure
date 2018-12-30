function [ n_timePlot ] = plotParameters( s_plot, s_protocol, s_current, ...
        s_constants, b_plotAll, n_delta, b_ignoreLastIndex )
%% PLOTPARAMETERS Plot the parameters found so far against the original values
%
% NOTE: The 'CurrentFigure' trick prevents the figure stealing the window focus.

    timer = tic;
    
    set(0, 'CurrentFigure', s_plot.figParam.fig);
    
    n_currentIdx = s_current.index;
    if (exist('b_ignoreLastIndex', 'var') && b_ignoreLastIndex)
        n_currentIdx = n_currentIdx - 1;
    else
        b_ignoreLastIndex = false;
    end
    if (b_plotAll)
        v_newData = 1 : n_currentIdx;
    else
        if (~ exist('n_delta', 'var'))
            n_delta = n_currentIdx - s_plot.lastPlotIndex;
        end
        v_newData = n_currentIdx - n_delta + 1 : n_currentIdx;
    end
    
    if (isempty(v_newData))
        % nothing to plot
        n_timePlot = toc(timer);
        return;
    end
    
    n_noParamsTotal = s_current.parametersCurrent.getTotalParamsNo();
    v_allParameters = [s_protocol(v_newData).parameters];
    m_allNewParams = zeros(n_noParamsTotal, size(v_allParameters, 2));
    for i = 1 : size(m_allNewParams, 2)
        m_allNewParams(:, i) = ...
            v_allParameters(1, i).getAllValuesIncludingTime();
    end
    v_allRefParams = ...
        s_constants.parametersReference.getAllValuesIncludingTime();
    n_plotIdx = 0;
    c_domains = s_constants.domains;
    for i = 1 : n_noParamsTotal
        switch (c_domains{i, 3})
            case {Constants.C_si_domainStandard, Constants.C_si_domainBoolean}
                n_plotIdx = n_plotIdx + 1;
                % nothing to do
            
            case Constants.C_si_domainConstant
                % skip constant parameters
                continue;
            
            otherwise
                error('Illegal domain');
        end
        subplotHandle = s_plot.figParam.subplots(n_plotIdx);
        set(s_plot.figParam.fig, 'CurrentAxes', subplotHandle);
        hold on
        
        % possible offset due to changing the designated switching parameter
        n_idx = s_plot.refMap(i);
        b_showReference = (n_idx ~= -1);
        if (b_showReference)
            y = v_allRefParams(n_idx);
        end
        
        % plot parameters
        children = get(subplotHandle, 'Children');
        if (isempty(children))
            % current parameters
            semilogy(v_newData, m_allNewParams(i, :), 'Color', 'b');
            
            % reference parameters
            if (b_showReference)
                semilogy([v_newData(1), v_newData(end)], [y y], 'r');
            end
        else
            % only update plot
            % NOTE: inverted order necessary
            
            if (b_showReference)
                n_childIdx = 2;
            else
                n_childIdx = 1;
            end
            
            % current parameters
            xdata = 1 : n_currentIdx;
            ydata = [get(children(n_childIdx), 'YData'), m_allNewParams(i, :)];
            set(children(n_childIdx), 'XData', xdata, 'YData', ydata);
            
            % reference parameters
            if (b_showReference)
                childXdata = get(children(1), 'XData');
                childXdata(2) = n_currentIdx;
                set(children(1), 'XData', childXdata);
            end
        end
        
        hold off
    end
    
    % update plots
    drawnow;
    
    n_timePlot = toc(timer);
end