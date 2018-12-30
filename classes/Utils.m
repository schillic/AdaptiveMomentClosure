classdef Utils
%UTILS Static helper class

    properties (Constant)
        % -- error handling --
        
        si_errorPrefix = 'PI: ';
    end
    
    methods (Static)
        function warn( si_msg )
        %% prints a tool-specific decorated warning message
        
            fprintf(sprintf('<strong>Warning: %s</strong>\n', si_msg));
        end
        
        function [ si_msg ] = wrapError( si_msg )
        %% adds a tool-specific decorator to an error message
        
            si_msg = [Utils.si_errorPrefix, si_msg];
        end
        
        function [ b_answer ] = isError( si_msg )
        %% adds a tool-specific decorator to an error message
        
            b_answer = (~ isempty(strfind(si_msg, Utils.si_errorPrefix)));
        end
        
        function odeSolverError( s_exception )
        %% handles ODE solver errors
        
            if (DEBUG.C_n_outputOdeExceptionAll > 0)
                if ((DEBUG.C_n_outputOdeExceptionAll == 2) && ...
                        Utils.isError(s_exception.message))
                    fprintf([s_exception.message, '\n']);
                else
                    fprintf(getReport(s_exception));
                end
            end
            if (DEBUG.C_n_haltOdeExceptionAll > 0)
                if ((DEBUG.C_n_haltOdeExceptionAll == 2) && ...
                        Utils.isError(s_exception.message))
                    % nothing to do
                else
                    keyboard; % intended
                end
            end
        end
    end
end