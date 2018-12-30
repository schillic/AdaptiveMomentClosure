function init( )
%INIT Initialization of the framework (call at the beginning of a session only)
% In order to initialize all things necessary for the program to work, this
% function needs to be called once at the beginning.

    % add StochDynTools toolbox
    addpath('functions', filesep, 'stochdyntools');
    
    % add Sundials toolbox
    sundialsPath = ['functions', filesep, 'sundialsTB'];
    addpath(sundialsPath);
    startup_STB(sundialsPath);
    
    % add examples subfolders
    addpath(genpath('examples'));
    
    % add functions and classes
    addpath('functions');
    addpath(['functions', filesep, 'inferenceCore']);
    addpath(['functions', filesep, 'simulation']);
    addpath(['functions', filesep, 'odeSolving']);
    addpath(['functions', filesep, 'utils']);
    addpath(['functions', filesep, 'plotting']);
    addpath(['functions', filesep, 'benchmarks']);
    addpath('classes');
    
    % initialize debug interruption
    debugInterrupt(1);
end