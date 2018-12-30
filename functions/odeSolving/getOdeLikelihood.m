function [ n_likelihood, s_odeMoments ] = getOdeLikelihood( ...
        s_current, s_constants, s_mdyn, si_closureName, s_simulationMoments )
%% GETODELIKELIHOOD Computes the ODE data and its likelihood wrt. the given data

    % get moments
    s_odeMoments = ...
        getOdeMoments(s_current, s_constants, s_mdyn, si_closureName);
    
    % compute likelihood
    if (isempty(s_simulationMoments))
        % simulation was skipped: do not assign a big value here, otherwise this
        % will always be the best parameter point
        n_likelihood = -Inf;
    else
        n_likelihood = ...
            getLikelihood(s_constants, s_odeMoments, s_simulationMoments);
        
        % only take the real part of the likelihood
        % In theory, the likelihood cannot be complex. Due to rounding errors,
        % however, this can happen when the determinant of the covariance matrix
        % becomes negative and hence its logarithm is not real anymore.
        if (~ isreal(n_likelihood))
            if (DEBUG.C_b_outputComplexLikelihood)
                fprintf(...
                    ' <strong>NOTE: likelihood is complex: %s</strong>\n', ...
                    num2str(n_likelihood, '%7.4f'));
                n_likelihood = real(n_likelihood);
            end
        elseif (isnan(n_likelihood))
            error(Utils.wrapError('Likelihood is NaN'));
        end
    end
end


function [ s_odeMoments ] = getOdeMoments( ...
        s_current, s_constants, s_mdyn, si_closureName )
%% extracts the moments from the ODE solver

    o_params = s_current.parametersCurrent;
    n_maxMeasurementIndex = length(s_constants.momentTimes);
    n_noTimeIntervals = o_params.getIntervalNo();
    if (n_noTimeIntervals == 1)
        s_odeMoments = [];
        b_overwrite = true;
    else
        m_odeMomentsInit = zeros(n_maxMeasurementIndex, ...
            length(s_constants.simulationIndices));
        
        s_odeMoments = struct(...
            'mean', m_odeMomentsInit, ...
            'variance', m_odeMomentsInit);
        if (s_constants.likelihoodMode == Constants.C_n_likelihoodTwo)
            % two measured species
            s_odeMoments.covariance = zeros(n_maxMeasurementIndex, 1);
        end
        b_overwrite = false;
    end
    v_x0 = s_mdyn.Mu.x0;
    n_measurementIdx = 1;
    [m_params, v_timePoints] = o_params.getEventValueMatrix();
    % for each time interval run the ODE solver
    for i = 1 : n_noTimeIntervals
        % measurement times
        [v_measurementTimesOde, n_startTime, n_ignoreEntries, ...
            c_timeFunctions] = getMeasurementTimes(...
            s_constants, v_timePoints, n_measurementIdx, i, b_overwrite);
        
        % parameter values
        v_params = m_params(:, i)';
        
        % TODO(T) is it sometimes reasonable to skip the iteration?
        % E.g., when there is more than one time parameter and they are almost
        % equal, like in [0, t1] [t1, t2], [t2, tmax] and the second interval is
        % (almost) empty?
        
        % run ODE solver in the specified order as long as one succeeds
        for j = 1 : length(s_constants.odeFunctions)
            n_odeSolverUsed = j;
            try
                % run ODE solver
                fun_odeSolver = s_constants.odeFunctions{j};
                m_odeData = fun_odeSolver(si_closureName, v_params, v_x0, ...
                    v_measurementTimesOde, n_startTime, c_timeFunctions, ...
                    s_constants.odeTerminationBounds(j));
                
                % no exception, accept this result
                break;
            catch exception
                if (DEBUG.C_b_outputOdeExceptionAny)
                    fprintf(getReport(exception));
                end
                if (DEBUG.C_b_haltOdeExceptionAny)
                    keyboard; % intended
                end
                % if there is no ODE solver left, throw the exception again
                if (j == length(s_constants.odeFunctions))
                    throw(exception);
                end
            end
        end
        
        % read moments from data
        [s_odeMomentsRaw, v_x0] = odeData2moments(s_constants, m_odeData, ...
            s_mdyn, n_odeSolverUsed, n_ignoreEntries);
        
        % collect moments for all time intervals
        [s_odeMoments, n_measurementIdx] = ...
            concatenateMoments(s_constants, s_odeMoments, s_odeMomentsRaw, ...
            n_measurementIdx, b_overwrite);
        
        % Early termination: All remaining time intervals are (almost) empty.
        % This effectively means the time intervals start after the
        % second to last measurement point. We can safely ignore the data for
        % such small intervals.
        if (n_measurementIdx > n_maxMeasurementIndex)
            break;
        end
    end
end


function [ v_momentTimes, n_startTime, n_ignoreEntries, ...
    c_timeFunctionsMod ] = getMeasurementTimes( s_constants, v_timePoints, ...
        n_lowerMeasurementIdx, n_timePointIdx, b_overwrite)
%% computes the measurement times (especially for time switches)

    v_momentTimes = s_constants.momentTimes;
    
    if (n_timePointIdx == 1)
        % first interval
        n_timeMin = 0;
    else
        % further intervals
        n_timeMin = v_timePoints(n_timePointIdx - 1);
    end
    if (n_timePointIdx == length(v_timePoints) + 1)
        % last interval
        n_timeMax = v_momentTimes(end);
    else
        % further intervals
        n_timeMax = v_timePoints(n_timePointIdx);
        assert(n_timeMax <= v_momentTimes(end), ...
            'There must not be any events after the last measurement time.');
    end
    
    if (~ b_overwrite)
        % compute upper measurement index
        if (n_timeMax == v_momentTimes(end))
            n_upperMeasurementIdx = length(v_momentTimes);
        else
            for n_upperMeasurementIdx = ...
                    n_lowerMeasurementIdx : length(v_momentTimes)
                if (v_momentTimes(n_upperMeasurementIdx) >= n_timeMax)
                    break;
                end
            end
        end
        
        % extract relevant measurement time frame
        v_momentTimes = s_constants.momentTimes(...
            n_lowerMeasurementIdx : n_upperMeasurementIdx);
    end
    
    % start time
    if (n_timePointIdx == 1)
        % the first measurement time is not necessarily the start time
        n_startTime = 0;
    else
        n_startTime = n_timeMin;
    end
    
    % determine whether some indices should be ignored
    if (isempty(v_momentTimes))
        % empty measurements: we take [point between start and end, end]
        v_momentTimes = ...
            [n_timeMax - (n_timeMax - n_startTime) / 2, n_timeMax];
        n_ignoreEntries = Constants.C_n_ignoreAll;
    elseif (size(v_momentTimes, 2) == 1)
        % The ODE solvers expect at least two measurement points,
        % otherwise they apply their own steps in between.
        % We cannot take the start time, as the first measurement must have
        % some distance to it. So we take the point between start and end.
        v_momentTimes = ...
            [v_momentTimes - (v_momentTimes - n_startTime) / 2, ...
            v_momentTimes];
        n_ignoreEntries = Constants.C_n_ignoreFirst;
    else
        n_ignoreEntries = Constants.C_n_ignoreNone;
    end
    
    % convert time functions to lambda functions of start time
    c_timeFunctionsMod = s_constants.timeFunctions;
    for i = length(c_timeFunctionsMod) : -1 : 1
        c_timeFunctionsMod{i} = @(t2)c_timeFunctionsMod{i}(n_startTime, t2);
    end
end


function [ s_moments, v_x0 ] = odeData2moments( ...
        s_constants, m_allData, s_mdyn, n_odeSolverUsed, n_ignoreEntries )
%% extracts mean/variance/covariance from the ODE data
% NOTE: We take the real part in case the numbers are complex.

    % some solvers return the data in a transposed way
    if (s_constants.odeTransposeData(n_odeSolverUsed))
        m_allData = m_allData';
    end
    
    if (~ s_constants.odeStartsAtFirstIndex(n_odeSolverUsed))
        % project away the first row (= initial condition)
        m_allData(1, :) = [];
    end
    
    % initial state = last state from previous run
    v_x0 = m_allData(end, :)';
    
    switch (n_ignoreEntries)
        case {Constants.C_n_ignoreNone, Constants.C_n_ignoreAll}
            % nothing to do
        
        case Constants.C_n_ignoreLast
            % remove the last row
            m_allData(end, :) = [];
            
        case Constants.C_n_ignoreFirst
            % remove the first row
            m_allData(1, :) = [];
        
        otherwise
            error(Utils.wrapError(...
                sprintf('Unknown mode: %d.', n_ignoreEntries)));
    end
    
    if (n_ignoreEntries == Constants.C_n_ignoreAll)
        m_moment1storage = [];
        m_moment2centered = [];
    else
        for i = length(s_constants.meanIndices) : -1 : 1
            n_meanIdx = s_constants.meanIndices(i);
            n_varianceIdx = s_constants.varianceIndices(i);
            
            % first moment
            v_moment1 = real(m_allData(:, n_meanIdx));
            % first moment for storage
            m_moment1storage(:, i) = v_moment1;
            % uncentered second moment
            v_moment2 = real(m_allData(:, n_varianceIdx));
            % centered second moment
            m_moment2centered(:, i) = v_moment2 - v_moment1.^2;
        end
    end
    s_moments = struct(...
        'mean', m_moment1storage, ...
        'variance', m_moment2centered);
    
    if (s_constants.likelihoodMode == Constants.C_n_likelihoodTwo)
        % two measured species
        
        % covariance
        v_covariance = real(m_allData(:, s_constants.covarianceIndex));
        % centered covariance
        v_covarianceCentered = v_covariance - ...
            m_moment1storage(:, 1) .* m_moment1storage(:, 2);
        
        s_moments.covariance = v_covarianceCentered;
    end
    
    % additional assertion that the indices are correct
    c_names = s_mdyn.Mu.mon;
    for i = length(s_constants.meanIndices) : -1 : 1
        n_meanIdx = s_constants.meanIndices(i);
        n_varianceIdx = s_constants.varianceIndices(i);
        
        v_indices = s_mdyn.Mu.ndx(:, n_meanIdx);
        assert(v_indices(n_meanIdx) == 1, 'Error in mean index.');
        assert(v_indices(n_varianceIdx) == 2, 'Error in variance index.');
        
        % names
        si_nameMean{i} = char(c_names(n_meanIdx));
        si_nameVar = char(c_names(n_varianceIdx));
        assert(strcmp(si_nameVar, sprintf('%s^2', si_nameMean{i})), ...
            'Error in variance name.');
        
        % covariance
        if (s_constants.likelihoodMode == Constants.C_n_likelihoodTwo)
            n_idxCov = s_constants.covarianceIndex;
            assert(v_indices(n_idxCov) == 1, 'Error in covariance index.');
            if (i == 1)
                si_nameCov = char(c_names(n_idxCov));
                assert(strcmp(si_nameCov, ...
                    sprintf('%s*%s', si_nameMean{1}, si_nameMean{2})) || ...
                    strcmp(si_nameCov, ...
                    sprintf('%s*%s', si_nameMean{2}, si_nameMean{1})), ...
                    'Error in covariance name.');
            end
        end
    end
end


function [ s_odeMoments, n_measurementIdxNew ] = concatenateMoments( ...
        s_constants, s_odeMoments, s_odeMomentsRaw, n_measurementIdxOld, ...
        b_overwrite )
%% concatenates the moments of iterative 
    
    if (b_overwrite)
        s_odeMoments = s_odeMomentsRaw;
        n_measurementIdxNew = n_measurementIdxOld;
    else
        n_measurementIdxNew = ...
            n_measurementIdxOld + length(s_odeMomentsRaw.mean);
        v_vector = (n_measurementIdxOld : (n_measurementIdxNew - 1));
        for i = length(s_constants.meanIndices) : -1 : 1
            s_odeMoments.mean(i, v_vector) = s_odeMomentsRaw.mean(i);
            s_odeMoments.variance(i, v_vector) = s_odeMomentsRaw.variance(i);
        end
        if (s_constants.likelihoodMode == Constants.C_n_likelihoodTwo)
            s_odeMoments.covariance(v_vector) = s_odeMomentsRaw.covariance;
        end
    end
end


function [ n_likelihood ] = getLikelihood( ...
        s_constants, s_newMoments, s_dataMoments )
%% computes the likelihood

    if (isempty(s_constants.momentsFilterIndices))
        s_newMomentsFiltered = s_newMoments;
    else
        c_fields = fieldnames(s_newMoments);
        s_newMomentsFiltered = struct();
        for i = 1 : length(c_fields)
            si_field = c_fields{i};
            s_newMomentsFiltered.(si_field) = ...
                s_newMoments.(si_field)(s_constants.momentsFilterIndices);
        end
    end
    switch (s_constants.likelihoodMode)
        case Constants.C_n_likelihoodSingle
            n_likelihood = 0;
            for i = 1 : length(s_constants.simulationIndices)
                n_likelihood = n_likelihood + ...
                    getLikelihoodSingle(s_newMomentsFiltered, s_dataMoments, i);
            end
        
        case Constants.C_n_likelihoodTwo
            n_measuredSpecies = length(s_constants.meanIndices);
            if (n_measuredSpecies == 1)
                n_likelihood = ...
                    getLikelihoodSingle(s_newMomentsFiltered, s_dataMoments, 1);
            elseif (n_measuredSpecies == 2)
                n_likelihood = ...
                    getLikelihoodTwo(s_newMomentsFiltered, s_dataMoments);
            else
                error(Utils.wrapError(...
                    ['More than two species are not allowed in the chosen ', ...
                    'likelihood mode.']));
            end
        
        otherwise
            error(Utils.wrapError('Unknown likelihood mode.'));
    end
end


function [ n_likelihood ] = getLikelihoodSingle( ...
        s_newMoments, s_dataMoments, n_speciesIndex )
%% computes the likelihood for a single species independent of others

    n_likelihood = 0;
    n_simulations = s_dataMoments.simulations;
    n_frac3 = (n_simulations - 3) / (n_simulations - 1);
    n_fracAll = 1 / n_simulations;
    v_dataM1 = s_dataMoments.mean(:, n_speciesIndex);
    v_dataM2 = s_dataMoments.variance(:, n_speciesIndex);
    v_dataM3 = s_dataMoments.moment3(:, n_speciesIndex);
    v_dataM4 = s_dataMoments.moment4(:, n_speciesIndex);
    v_newM1 = s_newMoments.mean(:, n_speciesIndex);
    v_newM2 = s_newMoments.variance(:, n_speciesIndex);
    
    for i = 1 : length(v_newM1)
        % if the mean is zero, we have nothing to measure
        if (v_dataM1(i) == 0)
            if (DEBUG.C_b_outputMeanZero)
                fprintf('mean = 0, skipping\n');
            end
            continue;
        end
        
        v_x = [v_newM1(i); v_newM2(i)];
        v_mu = [v_dataM1(i); v_dataM2(i)];
        v_xMinusMu = v_x - v_mu;
        
        m_sigma = n_fracAll * ...
            [v_dataM2(i), v_dataM3(i); ...
            v_dataM3(i), v_dataM4(i) - (n_frac3 * v_dataM2(i)^2)];
        
        n_difference = ...
            (log(det(m_sigma)) + v_xMinusMu' * (m_sigma \ v_xMinusMu)) / 2;
        n_likelihood = n_likelihood - n_difference;
    end
end


function [ n_likelihood ] = getLikelihoodTwo( s_newMoments, s_dataMoments )
%% computes the likelihood for two species

    n_likelihood = 0;
    n_simulations = s_dataMoments.simulations;
    n_nMinusOne = (n_simulations - 1);
    n_frac1 = 1 / n_nMinusOne;
    n_frac2 = (n_simulations - 2) / n_nMinusOne;
    n_frac3 = (n_simulations - 3) / n_nMinusOne;
    n_fracAll = 1 / n_simulations;
    v_dataM1_x = s_dataMoments.mean(:, 1);
    v_dataM2_x = s_dataMoments.variance(:, 1);
    v_dataM3_x = s_dataMoments.moment3(:, 1);
    v_dataM4_x = s_dataMoments.moment4(:, 1);
    v_dataM1_y = s_dataMoments.mean(:, 2);
    v_dataM2_y = s_dataMoments.variance(:, 2);
    v_dataM3_y = s_dataMoments.moment3(:, 2);
    v_dataM4_y = s_dataMoments.moment4(:, 2);
    v_dataM2_xy = s_dataMoments.covariance;
    v_dataM3_x2y = s_dataMoments.moment3_x2y;
    v_dataM3_xy2 = s_dataMoments.moment3_xy2;
    v_dataM4_x3y = s_dataMoments.moment4_x3y;
    v_dataM4_x2y2 = s_dataMoments.moment4_x2y2;
    v_dataM4_xy3 = s_dataMoments.moment4_xy3;
    v_newM1_x = s_newMoments.mean(:, 1);
    v_newM2_x = s_newMoments.variance(:, 1);
    v_newM1_y = s_newMoments.mean(:, 2);
    v_newM2_y = s_newMoments.variance(:, 2);
    v_newM2_xy = s_newMoments.covariance;
    
    for i = 1 : length(v_dataM1_x)
        n_mu1x = v_dataM1_x(i);
        n_mu1y = v_dataM1_y(i);
        
        % if the mean is zero, we have nothing to measure
        if ((n_mu1x == 0) || (n_mu1y == 0))
            if (DEBUG.C_b_outputMeanZero)
                fprintf('mean = 0, skipping\n');
            end
            continue;
        end
        
        n_mu2x = v_dataM2_x(i);
        n_mu2y = v_dataM2_y(i);
        n_mu2xy = v_dataM2_xy(i);
        n_mu3x = v_dataM3_x(i);
        n_mu3y = v_dataM3_y(i);
        n_mu3x2y = v_dataM3_x2y(i);
        n_mu3xy2 = v_dataM3_xy2(i);
        n_mu4x = v_dataM4_x(i);
        n_mu4y = v_dataM4_y(i);
        n_mu4x3y = v_dataM4_x3y(i);
        n_mu4x2y2 = v_dataM4_x2y2(i);
        n_mu4xy3 = v_dataM4_xy3(i);
        
        v_x = [v_newM1_x(i); v_newM1_y(i); v_newM2_x(i); v_newM2_xy(i); ...
            v_newM2_y(i)];
        v_mu = [n_mu1x; n_mu1y; n_mu2x; n_mu2xy; n_mu2y];
        v_xMinusMu = v_x - v_mu;
        
        n_entry1 = n_mu4x - n_frac3 * n_mu2x^2;
        n_entry2 = n_mu4x3y - n_frac3 * n_mu2x * n_mu2xy;
        n_entry3 = n_mu4x2y2 - n_frac3 * n_mu2x * n_mu2y;
        n_entry4 = ...
            n_mu4x2y2 - n_frac2 * n_mu2xy^2 + n_frac1 * n_mu2x * n_mu2y;
        n_entry5 = n_mu4xy3 - n_frac3 * n_mu2xy * n_mu2y;
        n_entry6 = n_mu4y - n_frac3 * n_mu2y^2;
        
        m_sigma = n_fracAll * ...
            [n_mu2x, n_mu2xy, n_mu3x, n_mu3x2y, n_mu3xy2; ...
            n_mu2xy, n_mu2y, n_mu3x2y, n_mu3xy2, n_mu3y; ...
            n_mu3x, n_mu3x2y, n_entry1, n_entry2, n_entry3; ...
            n_mu3x2y, n_mu3xy2, n_entry2, n_entry4, n_entry5; ...
            n_mu3xy2, n_mu3y, n_entry3, n_entry5, n_entry6];
        
        n_difference = ...
            (log(det(m_sigma)) + v_xMinusMu' * (m_sigma \ v_xMinusMu)) / 2;
        n_likelihood = n_likelihood - n_difference;
    end
end