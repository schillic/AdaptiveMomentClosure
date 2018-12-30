function [ s_current, n_timeMcmc ] = newParameters( ...
        s_current, s_constants, n_timeMcmc )
%% NEWPARAMETERS Compute new parameter values from the previous parameter values

    timer_mcmc = tic;
    
    o_oldParameterVector = s_current.parametersAccepted;
    v_oldParameters = o_oldParameterVector.getAllValuesIncludingTime();
    v_logOldParameters = log(v_oldParameters);
    n_parameters = length(v_oldParameters);
    v_newParameters = zeros(n_parameters, 1);
    
    b_chooseAgain = true;
    while (b_chooseAgain)
        b_chooseAgain = false;
        n_forwardProbability = 1;
        n_backwardProbability = 1;
        v_booleans = [];
        
        for i = 1 : n_parameters
            c_domain = s_constants.domains(i, :);
            
            % check for domain violation
            if (~ isempty(c_domain))
                % special domains handling
                switch (c_domain{1, 3})
                    case Constants.C_si_domainStandard
                        % nothing to do; continue below
                    
                    case Constants.C_si_domainConstant
                        % constant: do not change this value
                        v_newParameters(i) = v_oldParameters(i);
                        continue;
                    
                    case Constants.C_si_domainBoolean
                        % Boolean domain: remember for later modification
                        v_newParameters(i) = v_oldParameters(i);
                        v_booleans(end + 1) = i; %#ok<AGROW>
                        continue;
                    
                    otherwise
                        error('Illegal domain type');
                end
            end
            
            n_newParameterVariance = s_constants.newParameterVariance;
            
            % standard MCMC step with lognormal distribution
            v_newParameters(i) = ...
                lognrnd(v_logOldParameters(i), n_newParameterVariance);
            
            % check for domain violation
            if (~ isempty(c_domain))
                n_min = c_domain{1, 1};
                n_max = c_domain{1, 2};
                
                if ((v_newParameters(i) < n_min) || ...
                        (v_newParameters(i) > n_max))
                    % domain violation
                    b_chooseAgain = true;
                    break;
                end
            end
            
            % compute forward and backward probability
            n_forwardProbability = n_forwardProbability * ...
                pdf('logn', v_newParameters(i), v_logOldParameters(i), ...
                    n_newParameterVariance);
            n_backwardProbability = n_backwardProbability * ...
                pdf('logn', v_oldParameters(i), log(v_newParameters(i)), ...
                    n_newParameterVariance);
        end
        
        if (b_chooseAgain)
            if (DEBUG.C_b_outputDomainViolation)
                fprintf('Parameters were out of domain, trying again...\n');
            end
            
            continue;
        end
        
        % extra handling for Boolean values
        % Policy: reassign exactly one entry with probability 50%
        % (avoids drastic changes in one step)
        if (~ isempty(v_booleans))
            n_idx = randi(length(v_booleans));
%             for n_idx = length(v_booleans) : -1 : 1
                v_newParameters(v_booleans(n_idx)) = (rand() > 0.5);
%             end
        end
    end
    
    s_current.parametersCurrent = ...
        ParameterVector(o_oldParameterVector, v_newParameters');
    s_current.forwardProb = n_forwardProbability;
    s_current.backwardProb = n_backwardProbability;
    n_timeDiff = toc(timer_mcmc);
    n_timeMcmc = n_timeMcmc + n_timeDiff;
end