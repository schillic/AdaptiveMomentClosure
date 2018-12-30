function [ m_data ] = momentsSundials( si_closureName, v_parameters, v_x0, ...
        v_measurementTimes, n_startTime, c_timeFunctions, n_timeThreshold )
%% MOMENTSSUNDIALS Sundials ODE solver

    s_options = CVodeSetOptions('RelTol', 1.e-10, 'LinearSolver', 'Dense');
    fun_closure = str2func(si_closureName);
    % hacky trick to pass parameters as single argument
    c_params = num2cell(v_parameters);
    
    % ODE function (lambda expression)
    if (isempty(c_timeFunctions))
        fun_closure = @(t, x)sundialsWrapper(fun_closure, x, c_params);
    else
        fun_closure = @(t, x)sundialsWrapperTime(fun_closure, x, c_params, ...
            c_timeFunctions, t);
    end
    
    % ODE solver call
    CVodeInit(fun_closure, 'Adams', 'Functional', n_startTime, v_x0, s_options);
    [status, v_tOut, m_data] = CVode(v_measurementTimes, 'Normal'); %#ok<ASGLU>
    
    % check that there are as many measurements available as specified
    if (size(v_tOut, 2) ~= size(v_measurementTimes, 2))
        error(Utils.wrapError(...
            ['The number of measurements from the ODE solver does not ', ...
            'agree with the specification.']));
    end
end
%% --- end of main function ---


%% --- helper functions ---

function [ y, b_flag, v__newData ] = sundialsWrapper( fun_closure, v_x, ...
        c_params )
%% intermediate helper for the Sundials API (needs two more outputs)

    y = fun_closure(v_x, c_params{:});
    b_flag = 0;
    v__newData = [];
end

function [ y, flag, new_data ] = sundialsWrapperTime( fun_closure, v_x, ...
        c_params, c_funcTime, n_t )
%% intermediate helper for the Sundials API with time function

    y = fun_closure(v_x, c_params{:}, c_funcTime{:}, n_t);
    flag = 0;
    new_data = [];
end