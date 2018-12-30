function evaluateClosures( si_name )
%% EVALUATECLOSURES Evaluates only the adaptive closure data
% For all runs of a given model the results are parsed and the closure
% statistics are printed.

    s_list = dir(['results', filesep, si_name, '*']);
    
    for i = 1 : length(s_list)
        s_element = s_list(i);
        if (~ s_element.isdir)
            continue;
        end
        si_fullName = s_element.name;
        fprintf('<strong>Closure evaluation of ''%s''</strong>', si_fullName);
        evaluateData(si_fullName, 0, 0, -1);
    end
end