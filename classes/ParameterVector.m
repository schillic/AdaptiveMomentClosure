classdef ParameterVector < matlab.mixin.Copyable
%PARAMETERVECTOR Parameter wrapper

    properties (SetAccess = private)
        c_params = {}; % cell of all parameters
        v_times; % map from event (positive integer) to time point (parameter)
        n_totalNo = 0; % total number of parameters (cached)
    end
    
    methods
        function [ this ] = ParameterVector( varargin )
        %% constructor
        
            if (nargin == 0)
                % default constructor, use standard values
            elseif (nargin == 1)
                arg = varargin{1};
                if (isa(arg, 'ParameterVector'))
                    % simple copy function thanks to superclass
                    this = copy(arg);
                else
                    error('Constructor: Unsupported input.');
                end
            elseif (nargin == 2)
                arg = varargin{1};
                if (isa(arg, 'ParameterVector'))
                    % copy constructor with old parameter vector for the
                    % structure and vector for the new values
                    
                    arg2 = varargin{2};
                    assert(isrow(arg2), ...
                        'The values must be given as a row vector.');
                    n_vecLength = length(arg2);
                    
                    % parameters
                    c_oldParams = arg.c_params;
                    n_oldParamsSize = length(c_oldParams);
                    c_newParams = cell(n_oldParamsSize, 1);
                    n_offset = 1;
                    for i = 1 : n_oldParamsSize
                        cs_oldParam = c_oldParams{i};
                        n_paramSize = cs_oldParam.getTotalNo();
                        v_values = arg2(n_offset : n_offset + n_paramSize - 1);
                        c_newParams{i} = cs_oldParam.copy(v_values);
                        n_offset = n_offset + n_paramSize;
                    end
                    this.c_params = c_newParams;
                    
                    % event time points
                    this.v_times = arg2(...
                        (n_vecLength - size(arg.v_times, 2) + 1) : n_vecLength);
                    
                    % number of parameters
                    this.n_totalNo = arg.n_totalNo;
                elseif (iscell(arg))
                    % initialize fields from input
                    this.c_params = arg;
                    
                    lc_n_totalNo = 0;
                    for i = 1 : length(arg)
                        lc_n_totalNo = lc_n_totalNo + arg{i}.getTotalNo();
                    end
                    
                    arg = varargin{2};
                    if (isrow(arg) || isempty(arg))
                        % initialize fields from input
                        this.v_times = arg;
                    else
                        error('Constructor: Unsupported input (2).');
                    end
                    lc_n_totalNo = lc_n_totalNo + size(this.v_times, 2);
                    
                    % make sure there is an entry for each event (otherwise set -1)
                    this.ensureEventSize();
                    
                    % cache total parameter number
                    this.n_totalNo = lc_n_totalNo;
                else
                    error('Constructor: Unsupported input (1).');
                end
            else
                error('Constructor: Unsupported input.');
            end
        end
        
        function [ n_times ] = getTimesNo( this )
        %% returns the time map
        
            n_times = size(this.v_times, 2);
        end
        
        function [ n_intervalsNo ] = getIntervalNo( this )
        %% returns the number of time intervals
        
            n_intervalsNo = this.getTimesNo() + 1;
        end
        
        function [ n_time ] = getTime( this, n_event )
        %% returns the time value given an event
        
            assert(isEvent(n_event), 'Illegal event.');
            n_time = this.v_times(1, n_event);
        end
        
        function [ m_matrix, v_timesSorted ] = getEventValueMatrix( this )
        %% matrix (parameter x event) whose (i,j)-entry is the parameter i value
        %% after event j; vector with sorted event times
        % For this we assume the current event times fixed.
        
            lc_c_params = this.c_params;
            m_backMap = this.getEventBackMapping();
            m_matrix = zeros(length(lc_c_params), size(this.v_times, 2) + 1);
            for i = 1 : length(lc_c_params)
                m_matrix(i, 1) = lc_c_params{i}.getValue();
            end
            for i = 1 : size(m_backMap, 1)
                n_event = m_backMap(i, 1);
                for j = 1 : length(lc_c_params)
                    n_value = lc_c_params{j}.getValue(n_event);
                    if (n_value == Constants.C_n_defaultNoEvent)
                        m_matrix(j, i + 1) = m_matrix(j, i);
                    else
                        m_matrix(j, i + 1) = n_value;
                    end
                end
            end
            
            v_timesSorted = m_backMap(:, 2)';
        end
        
        function [ n_totalNo ] = getTotalParamsNo( this )
        %% returns the total number of parameters, including switches and time
        
            n_totalNo = this.n_totalNo;
        end
        
        function [ v_allParams ] = getAllValuesIncludingTime( this )
        %% get a vector containing all possible parameters
        % Those are both the parameter values for all time intervals and all the
        % time points.
        %
        % NOTE: The order of the vector entries is important and must not be
        % changed.
        
            v_allParams = zeros(this.n_totalNo, 1);
            lc_c_params = this.c_params;
            idx = 1;
            for i = 1 : length(lc_c_params)
                v_new = lc_c_params{i}.getAllValues();
                newIndices = idx : (idx + length(v_new) - 1);
                v_allParams(newIndices) = v_new;
                idx = idx + length(v_new);
            end
            v_allParams(idx : end) = this.v_times;
        end
        
        function [ c_names ] = getNameVector( this, v_symNames )
        %% returns a mapping from parameter index to name
        
            lc_n_totalNo = this.getTotalParamsNo();
            c_names = cell(lc_n_totalNo, 1);
            
            % time names
            n_timesNo = this.getTimesNo();
            if (n_timesNo == 1)
                c_names{lc_n_totalNo} = 't';
            elseif (n_timesNo > 1)
                for i = n_timesNo : -1 : 1
                    c_names{lc_n_totalNo - i + 1} = ...
                        ['t_', num2str(n_timesNo - i)];
                end
            end
            
            % parameter names
            n_paramsNo = length(this.c_params);
            idx = 1;
            for i = 1 : n_paramsNo
                v_values = this.c_params{i}.getAllValues();
                nextIdx = idx + length(v_values);
                si_name = char(v_symNames(i));
                
                % single parameter
                if (length(v_values) == 1)
                    c_newNames = {si_name};
                else
                    c_newNames = cell(length(v_values), 1);
                    for j = 1 : length(v_values)
                        c_newNames{j} = [si_name, '_', num2str(j)];
                    end
                end
                c_names(idx : nextIdx - 1) = c_newNames;
                idx = nextIdx;
            end
        end
    end
    
    methods (Access = private)
        function ensureEventSize( this )
        %% makes sure there is an entry for each event
        
            n_noEvents = size(this.v_times, 2);
            for i = 1 : length(this.c_params)
                this.c_params{i}.setValuesNo(n_noEvents);
            end
        end
        
        function [ b_answer ] = isEvent( this, n_event )
        %% checks whether the input is a valid event
        
            b_answer = true;
            if (~ isinteger(n_event))
                b_answer = false;
            elseif ((0 < n_event) && (n_event <= size(this.v_times, 2)))
                b_answer = false;
            end
        end
        
        function [ m_backMap ] = getEventBackMapping( this )
        %% matrix which maps sorted events to original ones
        % matrix (event x 2) whose i-th row contains (e_j,t_j) where t_j is the
        % i-th event time in sorted order and e_j is the associated event index
        
            n_noEvents = size(this.v_times, 2);
            m_backMap = zeros(n_noEvents, 2);
            
            % fill map
            for i = 1 : n_noEvents
                m_backMap(i, :) = [i, this.v_times(1, i)];
            end
            
            % sort indices (simple quadratic method)
            for i = 1 : n_noEvents
                for j = (i + 1) : n_noEvents
                    if (m_backMap(i, 2) > m_backMap(j, 2))
                        v_tmp = m_backMap(i, :);
                        m_backMap(i, :) = m_backMap(j, :);
                        m_backMap(j, :) = v_tmp;
                    end
                end
            end
        end
    end
end