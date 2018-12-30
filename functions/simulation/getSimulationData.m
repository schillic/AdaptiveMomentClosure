function [ c_simulationData, n_timeTotal ] = getSimulationData( ...
        s_constants, o_parameters, n_simulations, n_timeTotal )
%% GETSIMULATIONDATA Runs simulations and tracks time

    if (DEBUG.C_b_outputSimInfo)
        fprintf('running simulation with N = %d...\n', n_simulations);
    end
    
    timer_simulation = tic;
    
    % run simulations
    c_simulationData = s_constants.simulationFunction(o_parameters, ...
        n_simulations, s_constants);
    
    n_timeDiff = toc(timer_simulation);
    n_timeTotal = n_timeTotal + n_timeDiff;
    
    if (DEBUG.C_b_outputSimInfo)
        fprintf('simulation took %.2f seconds\n', n_timeDiff);
    end
end