function benchmarkEpsilon( )
%% BENCHMARKEPSILON Runs examples with different values for epsilon

    % activate epsilon evaluation mode
    global G_b_epsilonEvaluation;
    G_b_epsilonEvaluation = true; %#ok<NASGU>
    
    % set values for epsilon
    global G_n_epsilonEvaluation;
    G_n_epsilonEvaluation = [10 5 2 1 0.5 0.3 0.2 0.1]; %#ok<NASGU>
    
    % run all benchmarks
    benchmarkAll();
    
    % deactivate mode again
    G_b_epsilonEvaluation = [];
    G_n_epsilonEvaluation = [];
end