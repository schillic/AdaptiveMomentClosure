function [ n_maxIndex ] = maxIndex( v_vector )
%% MAXINDEX Finds the index of the maximum entry in a vector
    
    n_maxIndex = 1;
    n_maxValue = v_vector(1);
    for i = 2 : length(v_vector)
        if (n_maxValue < v_vector(i))
            n_maxIndex = i;
            n_maxValue = v_vector(i);
        end
    end
end