function [ s_moments, n_timeTotal ] = getSimulationMoments( ...
        s_constants, o_parameters, n_simulations, n_timeTotal, c_data )
%% GETSIMULATIONMOMENTS Converts the simulation data to centered moments
% If this function is called with data, it computes the moments for the data.
% Otherwise it first creates new data by running simulations.
% 
% NOTE: We use slightly imprecise formulae.
% For instance, the variance is given by 1/(n-1) (sum X^2 - 1/n (sum X)^2).
% Instead, we use E(X^2) - E(X)^2 = 1/n (sum X^2 - 1/n (sum X)^2).
% The difference is the first factor 1/(n-1) vs. 1/n.
% However, in our application these numbers are almost equal for sufficiently
% large values of n.

    if (isempty(c_data))
        % run simulations
        [c_data, n_timeTotal] = getSimulationData(s_constants, o_parameters, ...
            n_simulations, n_timeTotal);
    end
    
    v_simulationIndices = s_constants.simulationIndices;
    for i = length(v_simulationIndices) : -1 : 1
        % read relevant species at relevant measurement points
        x = c_data{v_simulationIndices(i)};
        x_sq = x.^2;
        x_cu = x.^3;
        x_tt4 = x.^4;
        assert(size(x, 2) == length(s_constants.measurementTimes), ...
                'The simulation should only return the correct measurements.');
        
        % terms for first species
        % E[X]
        Ex = mean(x, 1)';
        % E[X^2]
        Ex2 = mean(x_sq, 1)';
        % E[X^3]
        Ex3 = mean(x_cu, 1)';
        % E[X^4]
        Ex4 = mean(x_tt4, 1)';
        % E[X]^2
        Ex_sq = Ex.^2;
        % E[X]^3
        Ex_cu = Ex.^3;
        % E[X]^4
        Ex_tt4 = Ex.^4;
        
        % centered moments for first species
        % E[X]
        m1X(:, i) = Ex;
        % E[(X - E[X])^2]
        % = E[X^2] - E[X]^2
        m2X(:, i) = Ex2 - Ex_sq;
        % E[(X - E[X])^3]
        % = E[X^3] - (3 * E[X]^2 * E[X]) + (2 * E[X]^3)
        m3X(:, i) = Ex3 - (3 * Ex2 .* Ex) + (2 * Ex_cu);
        % E[(X - E[X])^4]
        % = E[X^4] - (4 * E[X^3] * E[X]) + (6 * E[X^2] * E[X]^2) - (3 * E[X]^4)
        m4X(:, i) = Ex4 - (4 * Ex3 .* Ex) + (6 * Ex2 .* Ex_sq) - (3 * Ex_tt4);
    end
    
    % data structure
    s_moments = struct(...
        'mean', m1X, ...
        'variance', m2X, ...
        'moment3', m3X, ...
        'moment4', m4X, ...
        'simulations', n_simulations);
    
    if (s_constants.likelihoodMode == Constants.C_n_likelihoodTwo)
        % additional data for two measured species
        % 
        % NOTE: We compute some numbers again which were already computed above.
        % 
        % NOTE: We exploit the fact that the species considered last is the
        % first species to be measured, so the species X above coincides with
        % species X here.
        
        y = c_data{v_simulationIndices(2)};
        y_sq = y.^2;
        y_cu = y.^3;
        assert(size(y, 2) == length(s_constants.measurementTimes), ...
            'The simulation should only return the correct measurements.');
        
        % E[Y]
        Ey = mean(y, 1)';
        % E[Y^2]
        Ey2 = mean(y_sq, 1)';
        % E[Y^3]
        Ey3 = mean(y_cu, 1)';
        % E[Y]^2
        Ey_sq = Ey.^2;
        % E[Y]^3
        Ey_cu = Ey.^3;
        % E[X * Y]
        Exy = mean(x .* y, 1)';
        % E[X^2 * Y]
        Ex2y = mean(x_sq .* y, 1)';
        % E[X * Y^2]
        Exy2 = mean(x .* y_sq, 1)';
        % E[X^3 * Y]
        Ex3y = mean(x_cu .* y, 1)';
        % E[X^2 * Y^2]
        Ex2y2 = mean(x_sq .* y_sq, 1)';
        % E[X * Y^3]
        Exy3 = mean(x .* y_cu, 1)';
        
        % mixed centered moments
        % covariance: E[(X - E[X]) * (Y - E[Y])]
        % = E[X * Y] - (E[X] * E[Y])
        m2XY = Exy - (Ex .* Ey);
        % E[(X - E[X])^2 * (Y - E[Y])]
        % = E[X^2 * Y] - (2 * E[X * Y] * E[X]) - (E[X^2] * E[Y]) ...
        %  + (2 * E[X]^2 * E[Y])
        m3X2Y = Ex2y - (2 * Exy .* Ex) - (Ex2 .* Ey) + (2 * Ex_sq .* Ey);
        % E[(X - E[X]) * (Y - E[Y])^2]
        % = E[X * Y^2] - (2 * E[X * Y] * E[Y]) - (E[X] * E[Y^2]) ...
        %  + (2 * E[X] * E[Y]^2)
        m3XY2 = Exy2 - (2 * Exy .* Ey) - (Ex .* Ey2) + (2 * Ex .* Ey_sq);
        % E[(X - E[X])^3 * (Y - E[Y])]
        % = E[X^3 * Y] - (3 * E[X] * E[X^2 * Y]) + (3 * E[X]^2 * E[X * Y]) ...
        %  - (3 * E[X]^3 * E[Y]) - (E[X^3] * E[Y]) + (3 * E[X^2] * E[X] * E[Y])
        m4X3Y = Ex3y - (3 * Ex .* Ex2y) + (3 * Ex_sq .* Exy) - ...
            (3 * Ex_cu .* Ey) - (Ex3 .* Ey) + (3 * Ex2 .* Ex .* Ey);
        % E[(X - E[X])^2 * (Y - E[Y])^2]
        % = E[X^2 * Y^2] - (2 * E[Y] * E[X^2 * Y]) - (2 * E[X] * E[X * Y^2]) ...
        %  + (4 * E[X] * E[Y] * E[X * Y]) + (E[X^2] * E[Y]^2) ...
        %  + (E[X]^2 * E[Y^2]) - (3 * E[X]^2 * E[Y]^2)
        m4X2Y2 = Ex2y2 - (2 * Ey .* Ex2y) - (2 * Ex .* Exy2) + ...
            (4 * Ex .* Ey .* Exy) + (Ex2 .* Ey_sq) + ...
            (Ex_sq .* Ey2) - (3 * Ex_sq .* Ey_sq);
        % E[(X - E[X]) * (Y - E[Y])^3] (analogous to m4X3Y)
        m4XY3 = Exy3 - (3 * Ey .* Exy2) + (3 * Ey_sq .* Exy) - ...
            (3 * Ey_cu .* Ex) - (Ey3 .* Ex) + (3 * Ey2 .* Ey .* Ex);
        
        % update data structure
        s_moments.covariance = m2XY;
        s_moments.moment3_x2y = m3X2Y;
        s_moments.moment3_xy2 = m3XY2;
        s_moments.moment4_x3y = m4X3Y;
        s_moments.moment4_x2y2 = m4X2Y2;
        s_moments.moment4_xy3 = m4XY3;
    end
end