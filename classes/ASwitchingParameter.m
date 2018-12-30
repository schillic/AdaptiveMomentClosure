classdef ASwitchingParameter < AParameter
%ASWITCHINGPARAMETER Parameter with several values for different time intervals
% When an event occurs (i.e., a certain time has passed), the value of this
% type of parameter can switch.

    properties
        v_values = []; % values for all time intervals
        v_intervalLengths = []; % interval lengths
    end
    
    methods
        function [ n_value ] = getValue( this, n_interval )
        %% returns the parameter value for a given interval
            
            for i = 1 : length(this.v_intervalLengths)
                n_interval = n_interval - this.v_intervalLengths(i);
                
                if (n_interval <= 0)
                    n_value = this.v_values(i);
                    return;
                end
            end
        end
        
        function [ v_values ] = getAllValues( this )
        %% returns all parameter values
            
            v_values = this.v_values;
        end
        
        function [ n_allParams ] = getAllParamsNo( this )
        %% returns the number of all parameter values
        
            n_allParams = length(this.v_values);
        end
        
        function [ v_intervals ] = getIntervals( this )
        %% returns the time intervals
        
            v_intervals = this.v_intervalLengths;
        end
        
        % -- functions to overwrite --
        
        function [ cs_newParameter ] = ...
                removeTimePoints( this, v_timePointsNeeded )
        %% removes time points
        
            assert(false, 'This method must be overwritten.');
            
            v_intervalLengthsNew = 1; % some initial value, will be overwritten
            n_currentIntervalIdx = 1;
            n_currentLength = this.v_intervalLengths(n_currentIntervalIdx);
            n_sum = 1;
            for i = 1 : length(v_timePointsNeeded)
                if (i == n_currentLength)
                    v_intervalLengthsNew(n_currentIntervalIdx) = n_sum;
                    n_currentIntervalIdx = n_currentIntervalIdx + 1;
                    if (n_currentIntervalIdx <= length(this.v_intervalLengths))
                        n_currentLength = n_currentLength + ...
                            this.v_intervalLengths(n_currentIntervalIdx);
                        n_sum = 1;
                    end
                elseif (v_timePointsNeeded(i))
                    n_sum = n_sum + 1;
                end
            end
            
            % construct new parameter
            cs_newParameter = ...
                ASwitchingParameter(this.v_values, v_intervalLengthsNew);
        end
        
        function [ cs_newParameter ] = copy( this, v_values )
        %% pseudo-copy constructor with new values
        
            assert(false, 'This method must be overwritten.');
            
            cs_newParameter = ...
                SwitchingParameter(v_values, this.v_intervalLengths);
        end
    end
end