function [ m_data ] = momentsMatlab( si_closureName, v_parameters, v_x0, ...
        v_measurementTimes, n_startTime, c_timeFunctions, n_timeThreshold )
%% MOMENTSMATLAB Matlab ODE solver for the given parameters and closure method

    s_options = odeset('RelTol', 1e-10);
    % termination control
    if (n_timeThreshold > 0)
        time_start = tic;
        s_options = odeset(s_options, ...
            'Events', ...
                @(t,y)controlTermination(t, y, n_timeThreshold, time_start));
    end
    
    fun_closure = str2func(si_closureName);
    % add the start time index as the first entry
    v_measurementTimes = [n_startTime, v_measurementTimes];
    % hacky trick to pass parameters as single argument
    c_params = num2cell(v_parameters);
    
    % ODE function (lambda expression)
    if (isempty(c_timeFunctions))
        fun_wrapped = @(t, x)fun_closure(x, c_params{:});
    else
        fun_wrapped = ...
            @(t, x)fun_closure(x, c_params{:}, c_timeFunctions{:}, t);
    end
    
    % ODE solver call
    [tOut, m_data] = ode45(fun_wrapped, v_measurementTimes, v_x0, s_options);
    
    % check that there are as many measurements available as specified
    if (size(tOut, 1) ~= size(v_measurementTimes, 2))
        error(Utils.wrapError(...
            ['The number of measurements from the ODE solver does not ' ...
            'agree with the specification.']));
    end
end


function [ n_value, b_terminal, n_direction ] = controlTermination( ...
        n_t, v_data, n_stop, time_start ) %#ok<INUSL>
%% controls termination of an ODE solver (fixed absolute run time bound)

    n_value = max(n_stop - toc(time_start), 0);
    b_terminal = 1;
    n_direction = -1;
end