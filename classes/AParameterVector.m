classdef (Abstract) AParameterVector < handle
%APARAMETERVECTOR Abstract parameter vector class

    methods (Abstract)
        %% returns the parameters
        [ c_params ] = getParams( this );
        
        %% returns the number of normal parameters
        [ n_params ] = getParamsNo( this );
        
        %% returns the sum of switching values for all normal parameters
        [ n_params ] = getParamsSwitchNo( this );
        
        %% returns the number of time intervals
        [ n_intervals ] = getIntervalsNo( this );
        
        %% get a vector containing the parameter values for all time intervals
        % For each parameter p the vector contains the values for each time
        % interval where the values differ. If p_i has m_i values for n time
        % intervals (m <= n), then the result contains sum_i m_i entries.
        %
        % NOTE: The order of the vector entries is important and must not be
        % changed.
        [ v_allParams ] = getAllValues( this );
        
        %% get a vector containing all possible parameters
        % Those are both the parameter values for all time intervals and all the
        % time points.
        %
        % NOTE: The order of the vector entries is important and must not be
        % changed.
        [ v_allParamsFinal ] = getAllValuesIncludingTime( this );
        
        %% get vector containing the parameter values for a given time
        [ v_values ] = getIntervalValues( this, n_time );
        
        %% get the parameter value for a given time
        [ n_value ] = getValue( this, n_param, n_time );
        
        %% get the time points
        %
        % NOTE: The order of the vector entries is important and must not be
        % changed.
        [ v_timePoints ] = getTimePoints( this );
        
        %% get a vector containing the number of values per parameter, including time
        [ v_offsets ] = getParamNoVector( this );
        
        %% get the parameter values and min/max time for a given time interval (index)
        [ v_values, n_timeMin, n_timeMax ] = getValuesAndTimeBounds( ...
                this, n_timeIntervalIdx, s_time );
        
        %% ensure time points are consistent
        % This method is used for changing a switching parameter dynamically.
        % Afterward, some of the time points may not be necessary anymore. Here
        % we remove them. The corrected parameters (where we also replace the
        % designated parameter) and time intervals are returned.
        [ c_params, v_timePoints, v_offsets ] = ensureConsistentTimePoints( ...
            this, n_idxIgnore, cs_paramInstead, n_timePointsRef );
    end
end