classdef DEBUG
%DEBUG Debug options

    properties (Constant)
        % -- output numbers --
        
        % output iteration number and time every n iterations (0 = deactivated)
        C_n_outputTotalIt = 100;
        
        % plot simulation iterations every n iterations (0 = deactivated) 
        C_n_outputSimIt = 1000;
        
        % plot mean and variance every n iterations (inf = deactivated)
        C_n_plotMeanVar = inf;
        
        % plot parameters every n iterations (inf = deactivated)
        C_n_plotParams = inf;
        
        % plot mode for data
        C_si_dataPlotMode = Constants.C_si_plotCircle;
        
        % -- output flags --
        
        % output simulation statistics
        C_b_outputSimInfo = 1;
        
        % output acceptance status
        C_b_outputAccepted = 0;
        
        % output final run time
        C_b_outputTotalTime = 1;
        
        % output result of switch updates (even for no updates)
        C_b_outputUpdateSwitchResult = 1;
        
        % output reevaluation information
        C_b_outputReevaluation = 1;
        
        % output if any ODE solver exception occurs
        C_b_outputOdeExceptionAny = 0;
        
        % output information when a closure crashed
        C_b_outputClosureCrash = 1;
        
        % output when likelihood is a complex number
        C_b_outputComplexLikelihood = 1;
        
        % output likelihood of all succeeding closures
        C_b_outputAllClosureLikelihoods = 1;
        
        % output information when choosing a closure
        C_b_outputClosureChoice = 1;
        
        % output when a parameter was proposed out of its domain
        C_b_outputDomainViolation = 1;
        
        % output when mean in reference data is zero
        C_b_outputReferenceMeanZero = 1;
        
        % output when mean is zero
        C_b_outputMeanZero = 0;
        
        % output when backtracking limit is reached
        C_b_outputBacktrackLimit = 1;
        
        % output when closure is loaded from file
        C_b_outputClosureLoad = 1;
        
        % output when data is generated or loaded from file
        C_b_outputDataFinished = 1;
        
        % output basic algorithm skeleton messages
        C_b_outputSkeleton = 1;
        
        % output results in the end
        C_b_outputResults = 1;
        
        % output when a parameter switch was updated
        C_b_outputUpdateSwitch = 1;
        
        % -- extended output flags --
        
        % output if all ODE solvers stop with an exception
        % 0 = deactive
        % 1 = activate
        % 2 = activate conditionally (only for non-intended errors)
        C_n_outputOdeExceptionAll = 2;
        
        % -- halt --
        
        % halt if debug figure is closed
        C_b_haltFigure = 1;
        
        % halt if all ODE solvers stop with an exception
        % 0 = deactive
        % 1 = activate
        % 2 = activate conditionally (only for non-intended errors)
        C_n_haltOdeExceptionAll = 2;
        
        % halt if any ODE solver exception occurs
        C_b_haltOdeExceptionAny = 0;
        
        % -- storage flags --
        
        % store results
        C_b_storeResults = 1;
        
        % -- extended storage flags --
        
        % write logs to a file
        % 0 = deactivate
        % 1 = activate and append
        % 2 = activate and overwrite
        C_n_storeLogs = 2;
        
        % -- plot flags --
        
        % show legend for plots
        C_b_showLegend = 0;
        
        % -- benchmark flags --
        
        % skip benchmarks which are flagged as working
        C_b_skipWorkingBenchmarks = 0;
    end
end