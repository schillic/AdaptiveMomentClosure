function [ s_closures, s_statistics, n_likelihood, s_odeMoments, ...
    b_crashed ] = runClosure( s_current, s_constants, s_closures, ...
        s_statistics, s_simulationMoments, n_closureIndex )
%% RUNCLOSURE Wrapper to run ODE solver for closure with exception handling

    timer_ode = tic;
    
    s_mdyn = s_closures.mdyns{n_closureIndex};
    si_closureName = s_constants.closureFileNames{n_closureIndex};
    b_crashed = false;
    
    try
        [n_likelihood, s_odeMoments] = getOdeLikelihood(s_current, ...
            s_constants, s_mdyn, si_closureName, s_simulationMoments);
    catch exception
        Utils.odeSolverError(exception);
        b_crashed = true;
        n_likelihood = -Inf;
        s_odeMoments = [];
        s_closures.crashes(end + 1) = n_closureIndex;
    end
    
    timeDiff = toc(timer_ode);
    
    % add time difference to the correct timer
    if (~ isempty(s_statistics))
        if (b_crashed)
            s_statistics.timeOdeCrashes = ...
                s_statistics.timeOdeCrashes + timeDiff;
        else
            s_statistics.timeOde = s_statistics.timeOde + timeDiff;
        end
    end
end