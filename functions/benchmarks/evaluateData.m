function evaluateData( ...
        si_folderName, b_plotParameters, b_plotMeanVariance, ...
        n_histogramIgnore, n_histogramRelDiff, b_storeFigures )
%% EVALUATEDATA Plotting of results written to files
% When a parameter search was run, the results are written to a new folder.
% These files are loaded and various evaluations are generated.

% check folder existence
si_folderPath = ['results', filesep, si_folderName, filesep];
if (~ exist(si_folderPath, 'dir'))
    error(Utils.wrapError(sprintf('Could not find folder ''%s''', ...
        si_folderPath)));
end

% variables to be loaded
c_fileNames = {'s_protocol', 's_current', 's_constants', 's_plot'};

% load data from files
for i = 1 : length(c_fileNames)
    try
        load([si_folderPath, c_fileNames{i}]);
    catch
        error(Utils.wrapError(...
            'Could not find file ''%s''', c_fileNames{i}));
    end
end

% table of parameter distance
fprintf('\n');
printParameterTable(s_current, s_constants);
fprintf('\n\n');

% closure statistics
printClosureTable(s_protocol, s_constants);
fprintf('\n');

% ask user for input if not specified
if (~ exist('b_plotParameters', 'var'))
    b_plotParameters = input(...
        'Plot parameters? [0/1]\n> ');
end
if (~ exist('b_plotMeanVariance', 'var'))
    b_plotMeanVariance = input(...
        'Plot mean and variance? [0/1]\n> ');
end
if (~ exist('n_histogramIgnore', 'var'))
    n_histogramIgnore = input(...
        ['How many indices should be ignored for the histogram? ', ...
        '[number, deactivate: -1]\n> ']);
end
b_plotHistogram = (n_histogramIgnore >= 0);
if (b_plotHistogram)
    if (~ exist('n_histogramRelDiff', 'var'))
        n_histogramRelDiff = input(...
            ['What is the relativ difference for the histogram on the X-axis? ', ...
            '[number > 0]\n> ']);
    end
end
if (~ exist('b_storeFigures', 'var'))
    b_storeFigures = true;
end

b_plot = (b_plotParameters || b_plotMeanVariance || b_plotHistogram);

% plot setup
if (b_plot)
    s_plotSetup = struct(...
        'parameters', b_plotParameters, ...
        'moments', b_plotMeanVariance, ...
        'momentsHidden', false, ...
        'histogram', b_plotHistogram ...
        );
    v_refMap = s_plot.refMap; %#ok<NODEF>
    s_plot = setupPlots(s_current, v_refMap, s_constants, s_plotSetup);
end

% plot parameters
if (b_plotParameters)
    plotParametersLocal(s_current, s_protocol, s_constants, s_plot, ...
        b_storeFigures);
    fprintf('\n');
end

% plot mean and variance for final parameters
if (b_plotMeanVariance)
    plotMeanVariance(s_plot, s_constants, s_constants.referenceMoments, ...
        Constants.C_n_plotData);
end

% compute histogram
if (b_plotHistogram)
    computeHistogram(s_protocol, s_constants, s_plot, n_histogramIgnore, ...    
        n_histogramRelDiff, b_storeFigures);
    fprintf('\n');
end

end
%% --- end of main function ---


%% --- helper functions ---

function printParameterTable( s_current, s_constants )
%% prints the table of parameter difference

    v_parametersReference = ...
        s_constants.parametersReference.getAllValuesIncludingTime();
    v_parametersBest = ...
        s_current.parametersBest.getAllValuesIncludingTime();
    
    if ((s_constants.parametersReference.getTimesNo() > 0) || ...
            (s_current.parametersBest.getTimesNo() > 0))
        error(Utils.wrapError('Time switch evaluation is not supported yet.'));
    end
    
    fprintf('parameter results:\nparam\t');
    for i = 1 : length(s_constants.parametersSymbolic)
        fprintf('%s\t', char(s_constants.parametersSymbolic(i)));
    end
    
    % reference parameters
    fprintf('\nV_true\t');
    for i = 1 : length(v_parametersReference)
        reference = v_parametersReference(i);
        fprintf('%s\t', printNumber(reference, 2));
    end
    
    % best parameters
    fprintf('\nV_best\t');
    for i = 1 : length(v_parametersBest)
        best = v_parametersBest(i);
        fprintf('%s\t', printNumber(best, 2));
    end
    
    % absolute difference
    fprintf('\nD_abs\t');
    for i = 1 : length(v_parametersReference)
        best = v_parametersBest(i);
        reference = v_parametersReference(i);
        difference = abs(best - reference);
        fprintf('%s\t', printNumber(difference, 2));
    end
    
    % relative difference
    fprintf('\nD_rel\t');
    for i = 1 : length(v_parametersReference)
        best = v_parametersBest(i);
        reference = v_parametersReference(i);
        difference = abs(best - reference);
        diffRelative = difference / reference * 100;
        fprintf('%s\t', printNumber(diffRelative, 2));
    end
end


function [ si_number ] = printNumber( n_number, n_floatingPointIndices )
%% prints a number in a nice way

    if (n_number > n_floatingPointIndices * 10^(- n_floatingPointIndices))
        si_number = sprintf('%.2f', n_number);
    else
        si_number = sprintf('%.2e', n_number);
    end
end


function printClosureTable( s_protocol, s_constants )
%% prints closure statistics

    c_closures = s_constants.closureMethods;
    
    if (size(c_closures, 1) == 1)
        fprintf('single closure (%s deg. %d) run with %d iterations\n', ...
            c_closures{1}, c_closures{2}, size(s_protocol, 1));
        b_singleClosure = true;
    else
        b_singleClosure = false;
    end
    
    s_closureDataRaw = getClosureDataRaw(s_protocol);
    try
        s_closureDataPerClosure = getClosureDataPerClosure(s_closureDataRaw);
    catch exception
        fprintf([exception.message, '\n']);
        return;
    end
    
    if (b_singleClosure)
        return;
    end
    
    fprintf('number of closure evaluations/switches: %d / %d\n', ...
        s_closureDataRaw.evaluations, s_closureDataRaw.switches);
    fprintf(['closure usage (abs/rel) for\n- each iteration (of %d) ', ...
        'and\n- for each evaluation (of %d):\nclosure\t'], ...
        size(s_protocol, 1), s_closureDataRaw.evaluations);
    for i = 1 : length(s_closureDataPerClosure)
        fprintf('%s%d\t', c_closures{i, 1}, c_closures{i, 2});
    end
    fprintf('total\nit_abs\t');
    for i = 1 : length(s_closureDataPerClosure)
        fprintf('%d\t', s_closureDataPerClosure(i).absolute);
    end
    fprintf('%d\nit_rel\t', size(s_protocol, 1));
    for i = 1 : length(s_closureDataPerClosure)
        fprintf('%.2f\t', s_closureDataPerClosure(i).percentage);
    end
    fprintf('---\nswi_abs\t');
    for i = 1 : length(s_closureDataPerClosure)
        fprintf('%d\t', s_closureDataPerClosure(i).choices);
    end
    fprintf('%d + 1 (init)\nswi_rel\t', s_closureDataRaw.switches);
    for i = 1 : length(s_closureDataPerClosure)
        fprintf('%.2f\t', s_closureDataPerClosure(i).choicesPercentage);
    end
    fprintf('---\n');
end


function [ s_closureDataRaw ] = getClosureDataRaw( s_protocol )
%% closure statistics helper 1

    n_last = 0;
    s_table = struct(...
        'closure', [], ...
        'length', []...
        );
    n_idx = 0;
    n_switches = 0;
    n_evaluations = 0;
    n_lastIdx = 1;
    i = 0;
    for n_current = [s_protocol.closure]
        i = i + 1;
        if (n_current ~= 0)
            n_evaluations = n_evaluations + 1;
            if (n_current ~= n_last)
                n_idx = n_idx + 1;
                s_table(n_idx).closure = n_current;
                if (n_idx > 1)
                    s_table(n_idx - 1).length = i - n_lastIdx;
                    n_switches = n_switches + 1;
                end
                n_last = n_current;
                n_lastIdx = i;
            end
        end
    end
    s_table(end).length = i - n_lastIdx + 1;
    s_closureDataRaw = struct(...
        'table', s_table, ...
        'total', size(s_protocol, 1), ...
        'evaluations', n_evaluations, ...
        'switches', n_switches...
        );
end


function [ s_closureDataPerClosure ] = getClosureDataPerClosure( ...
        s_closureDataRaw )
%% closure statistics helper 2

    s_table = s_closureDataRaw.table;
    
    % collect total amount of usage per closure
    s_closureDataPerClosure = struct(...
        'length', [], ...
        'choices', [], ...
        'absolute', [], ...
        'percentage', []...
        );
    if (isempty(s_table(1).closure))
        error('Parameter search had crashed.');
    end
    for i = 1 : length(s_table)
        n_closureIdx = s_table(i).closure;
        n_currentLength = s_table(i).length;
        if ((length(s_closureDataPerClosure) >= n_closureIdx) && ...
                (~ isempty(s_closureDataPerClosure(n_closureIdx).length)))
            s_closureDataPerClosure(n_closureIdx).length = ...
                s_closureDataPerClosure(n_closureIdx).length + n_currentLength;
            s_closureDataPerClosure(n_closureIdx).choices = ...
                s_closureDataPerClosure(n_closureIdx).choices + 1;
        else
            s_closureDataPerClosure(n_closureIdx).length = n_currentLength;
            s_closureDataPerClosure(n_closureIdx).choices = 1;
        end
    end
    
    n_perHundred = 100 / s_closureDataRaw.total;
    n_closureChoices = s_closureDataRaw.switches + 1;
    for i = 1 : length(s_closureDataPerClosure)
        if (~ isempty(s_closureDataPerClosure(i).length))
            s_closureDataPerClosure(i).absolute = ...
                s_closureDataPerClosure(i).length;
            s_closureDataPerClosure(i).percentage = ...
                n_perHundred * s_closureDataPerClosure(i).length;
            s_closureDataPerClosure(i).choicesPercentage = ...
                s_closureDataPerClosure(i).choices / n_closureChoices * 100;
        else
            s_closureDataPerClosure(i).absolute = 0;
            s_closureDataPerClosure(i).percentage = 0;
            s_closureDataPerClosure(i).length = 0;
            s_closureDataPerClosure(i).choices = 0;
            s_closureDataPerClosure(i).choicesPercentage = 0;
        end
    end
end


function plotParametersLocal( ...
        s_current, s_protocol, s_constants, s_plot, b_storeFigures )
%% plots parameters

    if (length(s_protocol) > 5000)
        si_longTimeNote = ', this may take a while';
    else
        si_longTimeNote = '';
    end
    fprintf('Plotting %d parameters%s...\n', ...
        length(s_protocol), si_longTimeNote);
    
    plotParameters(s_plot, s_protocol, s_current, s_constants, true);
    
    % save figure
    if (b_storeFigures)
        set(0, 'CurrentFigure', s_plot.figParam.fig);
        storeFigure(['results', filesep, 'parameter_search']);
    end
end


function computeHistogram( s_protocol, s_constants, s_plot, n_ignoreEntries, ...
        n_histogramRelDiff, b_storeFigures )
%% computes and plots a histogram

    v_allParameters = [s_protocol(:).parameters];
    assert(~ isempty(v_allParameters), ...
        'There should be at least one measurement.');
    n_parameters = v_allParameters(1).getTotalParamsNo();
    
    % ignore some entries
    if (n_ignoreEntries >= size(v_allParameters, 2))
        error(Utils.wrapError(...
            'There are not enough measurements to be ignored.'));
    end
    v_allParameters(1 : n_ignoreEntries) = [];
    
    n_measurements = size(v_allParameters, 2);
    
    if (n_measurements > 5000)
        si_longTimeNote = ', this may take a while';
    else
        si_longTimeNote = '';
    end
    fprintf(...
        'Reading %d parameters for histogram%s...\n', ...
        n_measurements, si_longTimeNote);
    
    m_allNewParams = zeros(n_parameters, n_measurements);
    for i = 1 : n_measurements
        if (mod(i, 1000) == 0)
            fprintf('%d, ', i);
            if (mod(i, 10000) == 0)
                fprintf('\n');
            end
        end
        m_allNewParams(:, i) = ...
            v_allParameters(1, i).getAllValuesIncludingTime();
    end
    if (n_measurements >= 1000)
        fprintf('\n');
    end
    
    s_figure = s_plot.figHist;
    set(0, 'CurrentFigure', s_figure.fig);
    v_parametersReference = ...
        s_constants.parametersReference.getAllValuesIncludingTime();
    for i = 1 : n_parameters
        subplotHandle = s_figure.subplots(i);
        set(s_figure.fig, 'CurrentAxes', subplotHandle);
        hold on
        
        % compute histogram
%         n_bins = 100;
%         hist(m_allNewParams(i, :), n_bins);
        [v_histogram, n_bins] = ComputeHistogram(m_allNewParams(i, :));
        
        % read reference parameter
        n_reference = v_parametersReference(i);
%         line([n_reference, n_reference], ...
%             get(subplotHandle, 'YLim'), ...
%             'Color', 'r');
        
        % relative parameters
        bar(n_bins/n_reference, v_histogram, 1);
        xlim([(1 - n_histogramRelDiff), (1 + n_histogramRelDiff)])
        ylim([0, max(v_histogram)])
    end
    
    % save figure
    if (b_storeFigures)
        set(0, 'CurrentFigure', s_plot.figHist.fig);
        storeFigure(['results', filesep, 'histogram']);
    end
end


function [v_histogram, v_bins] = ComputeHistogram(v_values, n_bins)
% computes histogram, i.e., the a posteriori distribution in bins

    if (nargin < 2)
        % sample standard deviation (of the a posteriori distribution)
        n_sigma = std(v_values);
        
        % optimal bin width (Scott's choice, see Wikipedia: Histogram)
        n_dOpt = 3.5 * n_sigma / (length(v_values))^(1/3);
        
        n_bins = (max(v_values) - min(v_values)) / n_dOpt;
    end
    
    % compute histogram;
    [v_histogram, v_bins] = hist(v_values, n_bins);
    n_deltaT = (max(v_bins) - min(v_bins)) / n_bins;
    v_histogram = v_histogram / sum(v_histogram * n_deltaT);
end