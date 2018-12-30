function benchmarkClosure( fun_example )
%% BENCHMARKCLOSURE Runs examples with different closures

    % activate single closure evaluation mode
    global G_b_singleClosureEvaluation;
    G_b_singleClosureEvaluation = true; %#ok<NASGU>
    
    if (exist('fun_example', 'var'))
        if (isa(fun_example, 'function_handle'))
            % run input function handle
            fun_example();
        elseif (ischar(fun_example) && exist(fun_example) ~= 5)
                % run input string interpreted as function
                fun_exec = str2func(fun_example);
                fun_exec();
        else
            error(Utils.wrapError('The input argument must be a function.'));
        end
    else
        % run all benchmarks
        benchmarkAll();
    end
    
    % deactivate mode again
    G_b_singleClosureEvaluation = false;
end