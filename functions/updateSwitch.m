function [ s_current, s_protocol, v_refMap ] = updateSwitch( ...
        s_current, s_protocol, n_switchParamIdx, n_iterationSwitch, n_seed )
%% UPDATESWITCH Updates the number of switches for the designated parameter

    % reset RNG
    resetRng(n_seed);

    if (n_switchParamIdx == 0)
        % update feature disabled, nothing to do
        v_refMap = (1 : s_current.parametersCurrent.getTotalParamsNo())';
        return;
    end
    
    error('The update feature is currently deactivated.');
    
    o_parameterVector = s_current.parametersCurrent;
    c_params = o_parameterVector.getParams();
    o_designatedParam = c_params{n_switchParamIdx};
    
    assert(isa(o_designatedParam, 'SwitchingParameter'), ...
        'Unsupported parameter type.');
    n_designatedSwitchesOld = o_designatedParam.getAllParamsNo() - 1;
    if (n_iterationSwitch == 0)
        % no switch; turn into a constant parameter
        o_designatedParam = ...
            ConstantParameter(o_designatedParam.getValue(1));
    else
        if (n_iterationSwitch == n_designatedSwitchesOld)
            % no change needed
        else
            % construct a new switching parameter
            
            % TODO(T) At the moment, only downgrading is supported
            assert((n_iterationSwitch < n_designatedSwitchesOld), ...
                'Upgrading the switch count is not supported yet.');
            
            n_remainingIndices = n_iterationSwitch + 1;
            v_values = o_designatedParam.getAllValues();
            v_intervals = o_designatedParam.getIntervals();
            n_sum = 0;
            for i = n_remainingIndices : length(v_intervals)
                n_sum = n_sum + v_intervals(i);
            end
            v_values = v_values(1 : n_remainingIndices);
            v_intervals = v_intervals(1 : n_remainingIndices);
            v_intervals(n_remainingIndices) = n_sum;
            
            o_designatedParam = SwitchingParameter(v_values, v_intervals);
        end
    end
    
    % check that the time points are still of the correct dimension
    c_paramsRef = s_current.parametersReference.getParams();
    n_designatedSwitchesRef = ...
        length(c_paramsRef{n_switchParamIdx}.getIntervals()) - 1;
    [c_params, v_timePoints, v_offsets] = ...
        o_parameterVector.ensureConsistentTimePoints(n_switchParamIdx, ...
        o_designatedParam, n_designatedSwitchesRef);
    
    % update data structures
    s_current.parametersCurrent = ParameterVector(c_params, v_timePoints);
    s_current.parametersAccepted = ParameterVector(c_params, v_timePoints);
    s_current.parametersBest = ParameterVector(c_params, v_timePoints);
    s_protocol(1).parameters = ParameterVector(c_params, v_timePoints);
    
    % mapping from new parameter indices to reference parameter indices
    n_newParams = s_current.parametersCurrent.getParamsSwitchNo();
    n_newTime = s_current.parametersCurrent.getIntervalsNo() - 1;
    n_totalNo = n_newParams + n_newTime;
    v_refMap = zeros(n_totalNo, 1);
    n_totalParamsNo = length(c_params);
    n_idx = 0;
    n_offset = 0;
    for i = 1 : n_totalParamsNo
        n_paramNo = c_params{i}.getAllParamsNo();
        for j = 1 : n_paramNo
            n_idx = n_idx + 1;
            v_refMap(n_idx) = n_idx + n_offset;
        end
        if (i == n_switchParamIdx)
            for j = n_designatedSwitchesOld : n_paramNo
                v_refMap(j) = -1;
            end    
            n_offset = -(n_paramNo - n_designatedSwitchesOld + 1) + ...
                getOffset(s_current, n_switchParamIdx, n_designatedSwitchesOld);
        end
    end
    % time parameters
    n_offsetOffset = v_refMap(n_newParams);
    for i = 1 : n_newTime
        n_newValue = v_offsets(i);
        if (n_newValue == -1)
            v_refMap(i + n_newParams) = -1;
        else
            v_refMap(i + n_newParams) = n_newValue + n_offsetOffset;
        end
    end
end
%% --- end of main function ---


%% --- helper functions ---


function [ n_offset ] = getOffset( s_current, n_switchParamIdx, n_designatedSwitchesOld )
%% calculate the offset

    c_paramsRef = s_current.parametersReference.getParams();
    n_offset = c_paramsRef{n_switchParamIdx}.getAllParamsNo() - ...
        n_designatedSwitchesOld + 1;
end