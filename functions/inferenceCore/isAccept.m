function [ b_answer ] = isAccept( n_likelihoodNew, s_current )
%% ISACCEPT Checks whether the parameters are accepted
% We compute the acceptance probability and compare to some random number.

    n_likelihoodOld = s_current.likelihoodPrevious;
    n_acceptProb = min(1, exp(n_likelihoodNew - n_likelihoodOld) * ...
                s_current.backwardProb / s_current.forwardProb);
    b_answer = (rand() < n_acceptProb);
end