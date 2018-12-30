classdef ConstantParameter < AParameter
%CONSTANTPARAMETER Parameter with one value for the whole time
% The value of this type of parameter does not change when an event occurs.

    properties (SetAccess = private)
        n_value = -1; % value of the parameter
    end
    
    methods
        function [ this ] = ConstantParameter( n_value )
        %% constructor
        
            if (nargin > 0)
                this.n_value = n_value;
            end
        end
        
        function [ n_value ] = getValue( this, unused )
        %% returns the parameter value (needs no further input for time)
        
            n_value = this.n_value;
        end
        
        function [ v_values ] = getAllValues( this )
        %% returns all parameter values
            
            v_values = this.n_value;
        end
        
        function [ n_allParams ] = getAllParamsNo( this )
        %% returns the number of all parameter values
        
            n_allParams = 1;
        end
        
        function [ v_intervals ] = getIntervals( this )
        %% returns the time intervals
        
            v_intervals = [];
        end
        
        function [ cs_newParameter ] = ...
                removeTimePoints( this, v_timePointsNeeded )
        %% removes time points
        
            % just return the old parameter
            cs_newParameter = this;
        end
        
        function [ cs_newParameter ] = copy( this, n_value )
        %% pseudo-copy constructor with new values
        
            cs_newParameter = ConstantParameter(n_value);
        end
    end
end