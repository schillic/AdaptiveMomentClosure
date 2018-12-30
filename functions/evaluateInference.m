function [ s_currentBest, s_protocolBest, b_continueSwitchInference ] = ...
    evaluateInference( s_currentNew, s_currentBest, s_protocolNew, ...
        s_protocolBest, f_likelihoodIncreaseFactor )
%% EVALUATEINFERENCE Evaluates the inference with switching parameters
    
    n_likelihoodBestNew = s_currentNew.likelihoodBest;
    fprintf('best likelihood: %.2f\n', n_likelihoodBestNew);
    
    if (isempty(s_currentBest))
        % no evaluation after first iteration necessary
        s_currentBest = s_currentNew;
        s_protocolBest = s_protocolNew;
        b_continueSwitchInference = true;
        return;
    end
    
    n_likelihoodBestOld = s_currentBest.likelihoodBest;
    
    % store data of newest iteration if the results were better
    if (n_likelihoodBestNew > n_likelihoodBestOld)
        s_currentBest = s_currentNew;
        s_protocolBest = s_protocolNew;
    end
    
    % check whether the next iteration should be tried
    % Multiply with specified improvement factor (special cases for negative
    % numbers).
    if (n_likelihoodBestOld >= 0)
        b_continueSwitchInference = (n_likelihoodBestNew > ...
            n_likelihoodBestOld * f_likelihoodIncreaseFactor);
    
    elseif (f_likelihoodIncreaseFactor > 0)
        b_continueSwitchInference = (n_likelihoodBestNew > ...
            n_likelihoodBestOld / f_likelihoodIncreaseFactor);
    
    else
        b_continueSwitchInference = (n_likelihoodBestNew > ...
            -n_likelihoodBestOld * f_likelihoodIncreaseFactor);
    end
        
    if (DEBUG.C_b_outputUpdateSwitchResult)
        if (n_likelihoodBestOld == 0)
            si_percent = '--';
        else
            si_percent = sprintf('%.2f', ...
                ((n_likelihoodBestNew - n_likelihoodBestOld) ...
                / abs(n_likelihoodBestOld) * 100));
        end
        fprintf('previously best likelihood: %.2f (improvement: %s %%)\n', ...
            n_likelihoodBestOld, si_percent);
    end
end