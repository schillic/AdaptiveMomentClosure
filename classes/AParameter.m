classdef (Abstract) AParameter < handle
%APARAMETER Abstract parameter class

    methods (Abstract)
        %% returns the value given a new event (-1 if no new value)
        %% If no event is passed, the initial value is returned.
        [ n_value ] = getValue( this, varargin );
        
        %% returns all parameter values
        [ v_values ] = getAllValues( this );
        
        %% returns the total number of parameter values
        [ n_totalNo ] = getTotalNo( this );
        
        %% pseudo-copy constructor with new values
        [ cs_newParameter ] = copy( this, v_values );
        
        %% makes sure there is an entry for each event (otherwise inserts -1)
        setValuesNo( this, n_newSize );
    end
end