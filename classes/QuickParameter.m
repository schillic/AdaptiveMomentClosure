classdef QuickParameter < AParameter
%% QUICKPARAMETER Quick and dirty parameter implementation
% Mapping from event to value (either -1 for no change or a non-negative value)

    properties
        n_initialValue = -1; % initial value
        v_values = []; % mapping from event (positive integer) to value
        v_valuesExtended = []; % extended mapping with entry for each event
    end
    
    methods
        function [ this ] = QuickParameter( varargin )
        %% constructor
        
            if (nargin == 0)
                % default constructor, use standard values
            elseif (nargin == 2)
                this.n_initialValue = varargin{1};
                this.v_values =  varargin{2};
            else
                error('Constructor: Unsupported input.');
            end
        end
        
        function [ n_value ] = getValue( this, varargin )
        %% returns the value given a new event (-1 if no new value)
        
            if (nargin == 1)
                n_event = 0;
            elseif (nargin == 2)
                n_event = varargin{1};
            else
                error('Illegal number of inputs.');
            end
            
            if (n_event == 0)
                n_value = this.n_initialValue;
            else
                assert(n_event <= length(this.v_valuesExtended), ...
                    'Illegal event.');
                n_value = this.v_valuesExtended(n_event);
            end
        end
        
        function [ v_values ] = getAllValues( this )
        %% returns all parameter values
            
            v_values = [this.n_initialValue, this.v_values];
        end
        
        function [ n_totalNo ] = getTotalNo( this )
        %% returns the total number of parameter values
        
            n_totalNo = 1 + length(this.v_values);
        end
        
        function [ cs_newParameter ] = copy( this, v_values )
        %% pseudo-copy constructor with new values
        
            cs_newParameter = QuickParameter(v_values(1), v_values(2 : end));
            cs_newParameter.setValuesNo(length(this.v_valuesExtended));
        end
        
        function setValuesNo( this, n_newSize )
        %% makes sure there is an entry for each event (otherwise inserts -1)
        
            n_oldSize = length(this.v_values);
            if (n_oldSize >= n_newSize)
                if (n_oldSize > n_newSize)
                    error('Constructor: Parameter has too many events.');
                end
                this.v_valuesExtended = this.v_values;
                return;
            end
            
            assert(n_oldSize < n_newSize, 'Can only add events.');
            for i = (n_newSize) : -1 : (n_oldSize + 1)
                this.v_valuesExtended(1, i) = Constants.C_n_defaultNoEvent;
            end
            if (n_oldSize > 0)
                this.v_valuesExtended(1, 1 : n_oldSize) = ...
                    this.v_values(1, 1 : n_oldSize);
            end
        end
    end
end