function [ b_answer ] = checkRecomputeClosure( s_current, s_constants, s_closures )
%% CHECKRECOMPUTECLOSURE Check whether the current closure must be reevaluated
% Typically, the closure must be reevaluated when the epsilon-ball was left.

    if (s_current.index == 1)
        % always run the evaluation the first time
        b_answer = true;
    elseif (s_closures.number == 1)
        % never reevaluate when only one closure remains
        b_answer = false;
    elseif (s_current.rejectionSeries == s_constants.rejectionLimit)
        % reevaluate when backtracking occurred too often
        if (DEBUG.C_b_outputBacktrackLimit)
            fprintf('backtracking limit reached, leaving loop\n');
        end
        b_answer = true;
    else
        % reevaluate when the relative epsilon ball was left, i.e., the
        % (maximum relative) distance to the center is too big
        
        v_center = s_current.parametersCenter.getAllValuesIncludingTime();
        v_current = s_current.parametersCurrent.getAllValuesIncludingTime();
        v_absDistance = v_center - v_current;
        n_relDistance = max(abs(v_absDistance ./ v_center));
        b_answer = (n_relDistance > s_constants.epsilon);
    end
end