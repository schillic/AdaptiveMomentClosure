function debugInterrupt( mode )
%% DEBUGINTERRUPT Called to initialize or check interruption condition
% argument mode:
% 1 = initialize
% 2 = check interruption condition and interrupt if it is satisfied
% 
% This was proposed here:
% http://stackoverflow.com/questions/3272541/ ...
% ... matlab-stop-and-continue-execution-from-debugger-possible/3273167#3273167

    global DEBUG_FIGURE;
    
    switch (mode)
        case 1
            % draw a small figure for debug interruption
            % Whenever the figure is closed and this condition is checked in the
            % code, the figure is recreated and execution is interrupted so the
            % user can inspect and change some variables. Note that it must be
            % continued explicitly, so the figure should not be closed
            % accidentally.
            DEBUG_FIGURE = getNewFigure();
            drawnow;
        case 2
            if (~ isempty(DEBUG_FIGURE))
                % wait a bit to check existence of figure
                pause(0.001);
                
                % if it does not exist anymore...
                if (~ ishandle(DEBUG_FIGURE))
                    % ... recreate it ...
                    DEBUG_FIGURE = getNewFigure();

                    % ... and interrupt
                    keyboard; % intended
                end
            end
    end
end


function [ newFigure ] = getNewFigure( )
%% creates a new debug figure

    newFigure = figure(100);
    set(newFigure, 'Position', [0, 0, 150, 100]);
end