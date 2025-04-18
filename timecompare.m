%% Compare computation time between data-driven and model-based approaches
function [X_data, X_data_model, time_model, time_data] = timecompare(AB1, AB2, Z0, U, W, N, Ad1, Bd1, Ad2, Bd2)


% Convert initial set (not timed)
Z0 = zono2zonotope(Z0);
Z0 = conZonotope(Z0);

% Set up parameters
totalsteps = N;

% Initialize cells for data-driven approach
X_data = cell(totalsteps+1, 1);
X_data{1} = Z0;

% Data-driven approach timing
total_time = 0;
for i = 1:totalsteps
    currentWidth = 2^(i - 1);
    for j = 1:currentWidth
        % Time reduction operation
        tic;
        X_data{i, j} = reduce(X_data{i, j}, 'girard', 400);
        total_time = total_time + toc;
        
        % Type conversion (not timed)
        X_data{i, j} = conZonotope2conZono(X_data{i, j});
        
        % Time halfspace intersection operations
        tic;
        X_data_negative = halfspaceIntersection(X_data{i, j}, [1 0], 0);
        X_data_positive = halfspaceIntersection(X_data{i, j}, [-1 0], 0);
        total_time = total_time + toc;
        
        % Type conversion (not timed)
        zonoNegZotope = conZono2conZonotope(X_data_negative);
        zonoPosZotope = conZono2conZonotope(X_data_positive);
        
        % Time cartesian product operations
        tic;
        cmz = cartProd(zonoNegZotope, U);
        cmz2 = cartProd(zonoPosZotope, U);
        total_time = total_time + toc;
        
        % Time matrix multiplication and addition operations
        tic;
        X_data{i + 1, 2 * j - 1} = MatzonotopetimesConzonotope(AB1, cmz) + W;
        X_data{i + 1, 2 * j} = MatzonotopetimesConzonotope(AB2, cmz2) + W;
        total_time = total_time + toc;
    end
end
time_data = total_time;

% Final reduction and conversion for data approach
for j = 1:2^(totalsteps)-1
    X_data{totalsteps+1, j} = reduce(X_data{totalsteps+1, j}, 'girard', 400);
    X_data{totalsteps+1, j} = conZonotope2conZono(X_data{totalsteps+1, j});
end

% Initialize cells for model-based approach
X_data_model = cell(totalsteps+1, 1);
X_data_model{1} = Z0;

% Model-based approach timing
total_time = 0;
for i = 1:totalsteps
    currentWidth = 2^(i - 1);
    for j = 1:currentWidth
        % Time reduction operation
        tic;
        X_data_model{i, j} = reduce(X_data_model{i, j}, 'girard', 400);
        total_time = total_time + toc;
        
        % Type conversion (not timed)
        X_data_model{i, j} = zonotope2zono(X_data_model{i, j});
        
        % Time halfspace intersection operations
        tic;
        X_data_negative = halfspaceIntersection(X_data_model{i, j}, [1 0], 0);
        X_data_positive = halfspaceIntersection(X_data_model{i, j}, [-1 0], 0);
        total_time = total_time + toc;
        
        % Type conversion (not timed)
        zonoNegZotope = zono2zonotope(X_data_negative);
        zonoPosZotope = zono2zonotope(X_data_positive);
        
        % Time matrix operations
        tic;
        X_data_model{i + 1, 2 * j - 1} = Ad1*zonoNegZotope + Bd1 * U + W;
        X_data_model{i + 1, 2 * j} = Ad2*zonoPosZotope + Bd2 * U + W;
        total_time = total_time + toc;
    end
end
time_model = total_time;

% Final conversion for model approach (not timed)
for j = 1:2^(totalsteps)
    X_data_model{totalsteps+1, j} = zonotope2zono(X_data_model{totalsteps+1, j});
end
end



