function [ s_plot ] = setupPlots( ...
        s_current, v_refMap, s_constants, s_plotSetup )
%% SETUPPLOTS Initialize plot figures if enabled

    if (s_constants.plot)
        o_parameters = s_current.parametersCurrent;
        n_totalValues = o_parameters.getTotalParamsNo();
        c_domains = s_constants.domains;
        for j = 1 : size(c_domains, 1)
            if (c_domains{j, 3} == Constants.C_si_domainConstant)
                n_totalValues = n_totalValues - 1;
            end
        end
        
        if (n_totalValues == 0)
            Utils.warn('There is no parameter to identify.');
        end
        
        % figure placement
        n_figurePlacement = s_constants.figurePlacement;
        n_figureIndex = 1;
        
        % parameter plot
        if (s_plotSetup.parameters)
            c_names = ...
                o_parameters.getNameVector(s_constants.parametersSymbolic);
            [fig_parameters, v_posP, n_figureIndex] = ...
                setupGeneralPlot(n_figureIndex);
            v_posLast = ...
                setupParameterPlot(fig_parameters, v_posP, 'parameters', ...
                n_figurePlacement);
            s_figParameters = setupSubplots(fig_parameters, c_names);
        else
            s_figParameters = [];
            v_posP = [1, 1, 0, 0];
            v_posLast = [];
        end
        
        % mean / variance / covariance plots
        if (s_plotSetup.moments || s_plotSetup.momentsHidden)
            n_measuredSpecies = length(s_constants.measuredVariableNames);
            c_names = cell(n_measuredSpecies * 2, 1);
            for i = 1 : n_measuredSpecies
                si_name = s_constants.measuredVariableNames{i};
                c_names{2 * i - 1} = sprintf('mean %s', si_name);
                c_names{2 * i} = sprintf('variance %s', si_name);
            end
            % covariance
            if (~ isempty(s_constants.covarianceIndex))
                c_names{end + 1} = sprintf('covariance %s / %s', ...
                    s_constants.measuredVariableNames{1}, ...
                    s_constants.measuredVariableNames{2});
            else
            end
            v_indices = [];
            if (s_plotSetup.moments)
                v_indices(end + 1) = 1;
            else
                s_figMeanVar = [];
            end
            if (s_plotSetup.momentsHidden)
                v_indices(end + 1) = 2;
            else
                s_figMeanVarRef = [];
            end
            for i = v_indices
                if (i == 1)
                    si_suffix = 'current';
                    si_visible = 'on';
                else
                    si_suffix = 'reference';
                    si_visible = 'off';
                end
                [fig_meanVar, v_posMV, n_figureIndex] = ...
                    setupGeneralPlot(n_figureIndex);
                v_posLast = setupMeanVariancePlot(fig_meanVar, v_posMV, ...
                    n_figurePlacement, v_posLast, [], si_suffix, si_visible);
                s_figure = setupSubplots(fig_meanVar, c_names);
                if (i == 1)
                    s_figMeanVar = s_figure;
                else
                    s_figMeanVarRef = s_figure;
                end
            end
        else
            s_figMeanVar = [];
            s_figMeanVarRef = [];
        end
        
        % histogram plot
        if (s_plotSetup.histogram)
            c_names = ...
                o_parameters.getNameVector(s_constants.parametersSymbolic);
            [fig_histogram, v_posH, n_figureIndex] = ...
                setupGeneralPlot(n_figureIndex); %#ok<NASGU>
            setupParameterPlot(fig_histogram, v_posH, 'histogram', ...
                n_figurePlacement, v_posP);
            s_figHist = setupSubplots(fig_histogram, c_names);
        else
            s_figHist = [];
        end
        
        % data structure
        s_plot = struct(...
            'isPlot', true, ...
            'lastPlotIndex', 0, ...
            'refMap', v_refMap, ...
            'figParam', s_figParameters, ...
            'figMeanVar', s_figMeanVar, ...
            'figMeanVarRef', s_figMeanVarRef, ...
            'figHist', s_figHist ...
            );
        
        % give focus to parameters plot
        figure(1);
    else
        s_plot = struct('isPlot', false);
    end
end
%% --- end of main function ---


%% --- helper functions ---


function [ n_rows, n_cols ] = getSquareGrid( n_totalValues )
%% constructs a near quadratic grid which contains all elements

    % find smallest square number which is greater or equal to total number
    n_rows = ceil(sqrt(n_totalValues));
    n_cols = n_rows;

    % shrink the grid by iteratively removing a row/column
    b_reduceRow = false;
    n_newSize = n_rows * n_cols;
    while (true)
        [b_reduceRow, n_difference] = ...
            findBigger(n_rows, n_cols, b_reduceRow);
        b_reduce = (n_newSize - n_difference >= n_totalValues);
        if (b_reduce)
            % reduce one row/column
            n_newSize = n_newSize - n_difference;
            if (b_reduceRow)
                n_rows = n_rows - 1;
            else
                n_cols = n_cols - 1;
            end
            b_reduceRow = false;
        else
            % undo last subtraction
            if (b_reduceRow)
                break;
            else
                b_reduceRow = true;
            end
        end
    end
end


function [ b_reduceRow, n_smaller ] = findBigger( n_rows, n_cols, b_preferRow )
%% find the bigger of rows and columns number
% The third argument is the tie breaker.

    if (b_preferRow)
        b_reduceRow = (n_rows <= n_cols);
    else
        b_reduceRow = (n_cols < n_rows);
    end
    
    if (b_reduceRow)
        n_smaller = n_cols;
    else
        n_smaller = n_rows;
    end
end


function [ v_pos ] = figurePlacementZero( v_pos )
%% placement of a figure in lower left corner of the screen

    v_pos = v_pos - [v_pos(1) - 1, v_pos(2) - 1, 0, 0];
end


function [ v_pos ] = figurePlacementRightOf( v_pos, v_posOther )
%% placement of a figure right of another one

    v_pos = [v_posOther(1) + v_posOther(3), ...
        v_posOther(2), ...
        v_pos(3), ...
        v_pos(4)];
end


function [ v_pos ] = figurePlacementAboveOf( v_pos, v_posOther )
%% placement of a figure above of another one

    v_pos = [v_posOther(1), ...
        v_posOther(2) + v_posOther(4) + Constants.C_n_figureBorderHeight, ...
        v_pos(3), ...
        v_pos(4)];
end


function [ fig_this, v_pos, n_figureIndex ] = setupGeneralPlot( n_figureIndex )
%% initializes a new figure

    fig_this = figure(n_figureIndex);
    clf(fig_this);
    v_pos = get(fig_this, 'OuterPosition');
    n_figureIndex = n_figureIndex + 1;
end

function [ v_posThis ] = setupParameterPlot( fig_this, v_posThis, si_name, ...
        n_figurePlacement, v_posLeft )
%% basic setup for parameter plot

    v_posThis = [v_posThis(1), v_posThis(2), ...
        Constants.C_n_plotParametersWidth, Constants.C_n_plotParametersHeight];
    switch (n_figurePlacement)
        case Constants.C_n_figureStandard
            % do nothing
        
        case Constants.C_n_figureZero
            v_posThis = figurePlacementZero(v_posThis);
        
        case Constants.C_n_figureGrid
            if (exist('v_posLeft', 'var'))
                v_posThis = figurePlacementRightOf(v_posThis, v_posLeft);
            else
                v_posThis = figurePlacementZero(v_posThis);
            end
        
        otherwise
            error(Utils.wrapError('Illegal placement setting.'));
    end
    set(fig_this, ...
        'Name', si_name, ...
        'OuterPosition', v_posThis);
end


function [ v_posThis ] = setupMeanVariancePlot( fig_this, ...
        v_posThis, n_figurePlacement, v_posLeft, v_posBelow, si_name, ...
        si_visible )
%% basic setup for mean/variance plot

    v_posThis = [v_posThis(1), v_posThis(2), ...
        Constants.C_n_plotWidth, Constants.C_n_plotHeight];
    switch (n_figurePlacement)
        case Constants.C_n_figureStandard
            % do nothing
        
        case Constants.C_n_figureZero
            v_posThis = figurePlacementZero(v_posThis);
        
        case Constants.C_n_figureGrid
            if (~ isempty(v_posLeft))
                v_posThis = figurePlacementRightOf(v_posThis, v_posLeft);
            end
            if (~ isempty(v_posBelow))
                v_posThis = figurePlacementAboveOf(v_posThis, v_posBelow);
            end
            if (isempty(v_posLeft) && isempty(v_posBelow))
                v_posThis = figurePlacementZero(v_posThis);
            end
        
        otherwise
            error(Utils.wrapError('Illegal placement setting.'));
    end
    
    set(fig_this, ...
        'Name', sprintf('moments %s', si_name), ...
        'Visible', si_visible, ...
        'OuterPosition', v_posThis);
end


function [ s_figure ] = setupSubplots( fig_this, c_names )
%% setup for subplots in a figure

    [n_rows, n_cols] = getSquareGrid(length(c_names));
    subplots = [];
    for i = length(c_names) : -1 : 1
        subplots(i) = subplot(n_rows, n_cols, i);
        set(fig_this, 'CurrentAxes', subplots(i));
        title(c_names{i});
    end
    s_figure = struct(...
        'fig', fig_this, ...
        'subplots', subplots);
end