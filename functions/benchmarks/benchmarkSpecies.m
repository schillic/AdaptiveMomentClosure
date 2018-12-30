function benchmarkSpecies( )
%% BENCHMARKSPECIES Runs examples with different species measurements

    % activate species measurement evaluation mode
    global G_b_speciesMeasuring;
    G_b_speciesMeasuring = true; %#ok<NASGU>
    
    % run all benchmarks
    benchmarkAll();
    
    % deactivate mode again
    G_b_speciesMeasuring = false;
end