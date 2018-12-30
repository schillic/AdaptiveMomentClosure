classdef SwitchingParameter < ASwitchingParameter
%SWITCHINGPARAMETER Parameter with several values for different time intervals
% When an event occurs (i.e., a certain time has passed), the value of this
% type of parameter can switch.
    
    methods
        function [ this ] = SwitchingParameter( v_values, v_intervalLengths )
        %% constructor
        
            if (nargin > 0)
                this.v_values = v_values;
                this.v_intervalLengths = v_intervalLengths;
                
                assert((length(this.v_values) == length(this.v_intervalLengths)), ...
                        'The vectors must have equal length.');
            end
        end
        
        function [ cs_newParameter ] = ...
                removeTimePoints( this, v_timePointsNeeded )
        %% removes time points
        
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
                SwitchingParameter(this.v_values, v_intervalLengthsNew);
        end
        
        function [ cs_newParameter ] = copy( this, v_values )
        %% pseudo-copy constructor with new values
        
            cs_newParameter = ...
                SwitchingParameter(v_values, this.v_intervalLengths);
        end
    end
end