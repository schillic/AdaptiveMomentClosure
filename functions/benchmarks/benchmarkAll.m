function benchmarkAll( )
%BENCHMARKALL Runs all examples

    if ((~ DEBUG.C_b_storeResults) && (~ DEBUG.C_n_storeLogs))
        error('Storing results is deactivated!');
    end
    if ((DEBUG.C_n_haltOdeExceptionAll ~= 0) || DEBUG.C_b_haltOdeExceptionAny)
        error('Halting must be deactivated!');
    end
    
    for i = 1 : 6
        switch (i)
            case 1
                fun_run = @runAphids;
                b_works = 1;
            
            case 2
                fun_run = @runJCP;
                b_works = 1;
            
            case 3
                fun_run = @runGeneTF;
                b_works = 1;
            
            case 4
                fun_run = @runPNAS;
                b_works = 1;
            
            case 5
                fun_run = @runJCP2;
                b_works = 0;
            
            case 6
                fun_run = @runTwoGenes;
                b_works = 1;
        end
        
        if (b_works && DEBUG.C_b_skipWorkingBenchmarks)
            continue;
        end
        
        clc;
        try
            fun_run();
        catch exception
            fprintf(getReport(exception));
        end
        pause(1);
    end
end