function [ s_constants, s_closures ] = fixClosureCrashes( ...
        s_constants, s_closures )
%% FIXCLOSURECRASHES Fix the result for closures which crashed
% There is a flag for removing the listed closures so they are never used again.

    v_crashIndices = s_closures.crashes;
    for i = length(v_crashIndices) : -1 : 1
        n_idxToRemove = v_crashIndices(i);
        
        if (DEBUG.C_b_outputClosureCrash)
            fprintf('<strong>closure %d (%s deg. %d) crashed', ...
                n_idxToRemove, s_closures.name{n_idxToRemove}, ...
                s_closures.degree(n_idxToRemove));
        end
        
        if (s_constants.removeCrashingClosures)
            fprintf(' and is disbanded</strong>\n');
            
            s_constants.closureFileNames(n_idxToRemove) = [];
            s_closures.name(n_idxToRemove) = [];
            s_closures.degree(n_idxToRemove) = [];
            s_closures.mdyns(n_idxToRemove) = [];
            
            % remove from likelihood vector (so it is not considered)
            if (~ isempty(s_closures.likelihoods))
                s_closures.likelihoods(n_idxToRemove) = [];
            end
        else
            fprintf('</strong>\n');
            
            % set likelihood to negative infinity (so it is not considered)
            if (~ isempty(s_closures.likelihoods))
                s_closures.likelihoods(n_idxToRemove) = -inf;
            end
        end
    end
    
    % offset to reference the correct plot index
    % NOTE: The index of crashed closures does not matter.
    v_idxOffset = zeros(length(s_closures.name), 1);
    if (s_constants.removeCrashingClosures)
        % reduce number of closures
        s_closures.num = s_closures.num - length(v_crashIndices);
    end
    n_offset = 0;
    n_next = v_crashIndices(1);
    for j = 1 : length(v_idxOffset)
        if (j == n_next)
            n_offset = n_offset + 1;
            if (n_offset < length(v_crashIndices))
                n_next = v_crashIndices(n_offset + 1);
            else
                n_next = 0;
            end
        end
        v_idxOffset(j) = n_offset;
    end
    s_closures.idxOffset = v_idxOffset;
end