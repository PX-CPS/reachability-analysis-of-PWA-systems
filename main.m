% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Example:
%   Data-driven reachability for a PWA system
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Parameters
clear;
close all;

N = 8;     % Number of discrete time steps
overapproximation_steps = 100; %overapproximation of family set steps
bounds = 3; % Bounds on state space for analysis (assuming [-bounds bounds] in each dimension)

% System dynamics
Ad1 = [0.75,0.25;-0.25,0.75];
Bd1= [-0.25;-0.25];
Ad2 = [0.75,-0.25;0.25,0.75];
Bd2 = [0.25;-0.25];

%% Compute the overapproximation of various model matrices
dim_x = size(Ad1,1);

X0 = zonotope(ones(dim_x,1)+ 1, 0.3 * diag(ones(dim_x,1)));
U_OverAPP = zonotope(-4, 0.025);

% Initialize results table for performance comparison
results = table('Size', [9, 3], 'VariableTypes', {'double', 'double', 'double'}, ...
                'VariableNames', {'N', 'Time_Data', 'Time_Model'});

% Find a good initial condition set
cx = [-0.255; 1.395];	% Crosses the guard once
Gx = 0.2*eye(2);
Z0 = zono(Gx,cx);
% Propagate backwards two steps
Ainv = Ad1^(-1);
Z0 = Ainv*(Z0 + -1*Bd1);
Z0 = Ainv*(Z0 + -1*Bd1);

% Create empty input set
U_new = conZonotope(1,0);

% overapproximation - family sets of truth model
W_noise = zonotope(zeros(dim_x,1), 0.0001*ones(dim_x,1));
[Ab1, Ab2] = OverApproximation(X0, U_OverAPP, W_noise, dim_x, Ad1, Bd1, Ad2, Bd2, overapproximation_steps);

% Time comparison for different numbers of steps
for i = 1:N+1
    [Z_OverApp, Z,time_model, time_data] = timecompare(Ab1, Ab2, Z0, U_new, W_noise, i, Ad1, Bd1, Ad2, Bd2); 
    % Store results
    results(i, :) = {i, time_data, time_model};
end

% Display results table
disp(results);

%% Visualization of performance comparison
% Create figure with specific size
figure('Position', [100, 100, 800, 600])

% Create the plot directly from results table
plot(results.N, results.Time_Data, 'b-o', 'LineWidth', 3, 'MarkerSize', 8, 'DisplayName', 'Computation Time From Data')
hold on
plot(results.N, results.Time_Model, 'r-x', 'LineWidth', 3, 'MarkerSize', 8, 'DisplayName', 'Computation Time From Model')

% Customize the plot
grid on
box on

% Set axis labels with LaTeX interpreter
xlabel('$Step$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Time (seconds)', 'Interpreter', 'latex', 'FontSize', 14)

% Customize the legend
legend('Location', 'northwest')
legend('boxoff')

% Set font size for axis numbers
set(gca, 'FontSize', 12)

% Add minor grid lines
grid minor

% Adjust axis limits with some padding
xlim([3.5 10.5])
ylim([0 max(results.Time_Data)*1.1])

% %% Exponential fitting for computation times
% % Prepare data
% N = results.N;
% Time_Data = results.Time_Data;
% Time_Model = results.Time_Model;
% 
% % Create exponential fit for both data sets
% ft = fittype('a*exp(b*x)', 'independent', 'x');
% [curve_data, gof_data] = fit(N, Time_Data, ft);
% [curve_model, gof_model] = fit(N, Time_Model, ft);
% 
% % Create denser point set for smooth curves
% N_fine = linspace(min(N), max(N), 100)';
% 
% % Create figure
% figure('Position', [100, 100, 800, 600])
% 
% % Plot original data points and fitted curves
% plot(N, Time_Data, 'b.', 'MarkerSize', 27, 'DisplayName', 'Data-Based')
% hold on
% plot(N, Time_Model, 'r.', 'MarkerSize', 27, 'DisplayName', 'Model-Based')
% plot(N_fine, curve_data(N_fine), 'b--', 'LineWidth', 3, 'DisplayName', sprintf('Data-Based Fit: %.5fe^{%.3fx}', curve_data.a, curve_data.b))
% plot(N_fine, curve_model(N_fine), 'r--', 'LineWidth', 3, 'DisplayName', sprintf('Model-Based Fit: %.5fe^{%.3fx}', curve_model.a, curve_model.b))
% 
% % Set figure properties
% grid on
% box on
% xlabel('N', 'FontSize', 14)
% ylabel('Time (seconds)', 'FontSize', 14)
% title('Exponential Fitting for Computation Time', 'FontSize', 16)
% legend('Location', 'northwest')
% 
% % Print fitting results
% fprintf('Over Approximation Fitting:\n');
% fprintf('f(x) = %.4f * e^(%.4f * x)\n', curve_data.a, curve_data.b);
% fprintf('R^2 = %.4f\n\n', gof_data.rsquare);
% 
% fprintf('Model Fitting:\n');
% fprintf('f(x) = %.4f * e^(%.4f * x)\n', curve_model.a, curve_model.b);
% fprintf('R^2 = %.4f\n', gof_model.rsquare);

%% Model-Based Approach
X0 = zono(bounds*eye(2),zeros(2,1));
X0a = halfspaceIntersection(X0,[1 0],0);
X0b = halfspaceIntersection(X0,[-1 0],0);

Phia = [eye(2);Ad1]*X0a + [zeros(2,1);Bd1];
Phib = [eye(2);Ad2]*X0b + [zeros(2,1);Bd2];

Phi = union(Phia,Phib);


figure; hold on;
title('Over Approximation of Various Conditions And True Model Results')
axis([-2 2 -1 3]);
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
set(gca,'fontsize',18,'fontname','times new roman')
p = plot([0 0], [-bounds bounds], 'g--','linewidth',2); % Guard

colors = interp1([1;N+1],[0 0 1;1 0 0],1:1:N+1);

plot(Z0,colors(1,:),1) % Plot initial condition set

% Compute and plot reachable sets
Z = Z0;
for i = 1:N
    Z = [zeros(2) eye(2)]*and(Phi,Z,[eye(2) zeros(2)]) + zonotope2zono(W_noise);
    plot(Z,colors(i,:),1) % Plot reachable set
end


%% Data-based Approach
% Z_OverApp = CalculatedDataFromOverApp(Ab1, Ab2, Z0, U_new, W_noise,steps);


plot(Z_OverApp{1,1},colors(N+1,:),0.1)
hold on % Plot reachable set

for i = 2:size(Z_OverApp,1) -1
    for j = 1:2^(i-1)
        if j <= 56 % less than 2^(i-1) to save time, plotting is so slow
            plot(Z_OverApp{i,j},colors(N+1,:),0.3) % Plot reachable set
        end
    end
end


legend(p,{'Guard'},'interpreter','latex')




% %% Compute next step sets from data
% function Data = CalculatedDataFromOverApp(AB1, AB2, Z0, U, W, steps)
% Z0 = zono2zonotope(Z0);
% Z0 = conZonotope(Z0);
% 
% % Set number of steps in analysis
% totalsteps = steps;
% X_data = cell(totalsteps+1,1);
% 
% % Initialize sets for loop
% X_data{1} = Z0;
% 
% for i=1:totalsteps
%     % Calculate sets for current step
%     currentWidth = 2^(i - 1);  % Current step width, number of j
%     for j = 1:currentWidth
%         % Convert X_data{i, j} to zonotope and reduce generators
%         X_data{i, j} = reduce(X_data{i, j}, 'girard', 400);
%         
%         % Convert to conZono type
%         X_data{i, j} = conZonotope2conZono(X_data{i, j});
%         
%         % Split set at x1 = 0
%         X_data_negative = halfspaceIntersection(X_data{i, j}, [1 0], 0);
%         X_data_positive = halfspaceIntersection(X_data{i, j}, [-1 0], 0);
%         
%         % Convert split parts back to zonotope
%         zonoNegZotope = conZono2conZonotope(X_data_negative);
%         zonoPosZotope = conZono2conZonotope(X_data_positive);
%         
%         % Cartesian product with input set
%         cmz = cartProd(zonoNegZotope, U);
%         cmz2 = cartProd(zonoPosZotope, U);
%         
%         % Convert matrix zonotopes to constraint matrix zonotopes
%         AB12 = matZonotope2conMatZonotope(AB1);
%         AB22 = matZonotope2conMatZonotope(AB2);
%         
%         % Calculate next states
%         X_data{i + 1, 2 * j - 1} = AB12*cmz + W;
%         if j == 8 && i== 4
%         
%         else
%           X_data{i + 1, 2 * j} = AB22*cmz2 + W;
%         end
% 
%     end
% end
% 
% % Convert final results to conZono
% for j = 1:2^(totalsteps)-1
%     X_data{totalsteps+1,j} = conZonotope2conZono(X_data{totalsteps+1,j});
% end
% 
% Data = X_data;
% end



%% Overapproximation- Family sets of truth model
function [AB1, AB2] = OverApproximation(X0, U, W, dim_x, A1, B1, A2, B2, overapproximation_steps)
% Function to compute over-approximations of system matrices

% Number of trajectories and time steps
initpoints = 1;
steps = overapproximation_steps;
totalsamples = initpoints * steps;

% Randomly choose constant inputs for each step
for i = 1:totalsamples
    u(i) = randPoint(U);
end

% Simulate the system to get the data
x0 = X0.center;
x(:, 1) = x0;
index = 1;
index_first = 1;
index_second = 1;
index_first_0 = 1; index_first_1 = 1;
index_second_0 = 1; index_second_1 = 1;

for j = 1:dim_x:initpoints * dim_x
    x(j:j + dim_x - 1, 1) = randPoint(X0);
    for i = 1:steps
        x1 = x(j, i);
        if (x1 <= 0)
            x_first(j:j + dim_x - 1, index_first) = x(j:j + dim_x - 1, i);
            x_meas_vec_first_0(:, index_first_0) = x_first(j:j+dim_x-1, index_first);
            x(j:j + dim_x - 1, i + 1) = A1 * x(j:j + dim_x - 1, i) + B1 * u(index) + randPoint(W);
            x_meas_vec_first_1(:, index_first_1) = x(j:j + dim_x - 1, i + 1);
            utraj_first(j, index_first) = u(index);
            index_first = index_first + 1;
            index_first_0 = index_first_0 + 1;
            index_first_1 = index_first_1 + 1;
        else
            x_second(j:j + dim_x - 1, index_second) = x(j:j + dim_x - 1, i);
            x_meas_vec_second_0(:, index_second_0) = x_second(j:j+dim_x-1, index_second);
            x(j:j + dim_x - 1, i + 1) = A2 * x(j:j + dim_x - 1, i) + B2 * u(index) + randPoint(W);
            x_meas_vec_second_1(:, index_second_1) = x(j:j + dim_x - 1, i + 1);
            utraj_second(j, index_second) = u(index);
            index_second = index_second + 1;
            index_second_0 = index_second_0 + 1;
            index_second_1 = index_second_1 + 1;
        end
        utraj(j, i) = u(index);
        index = index + 1;
    end
end

% Organize data for each subsystem
X_first_0T = x_meas_vec_first_0;
X_first_1T = x_meas_vec_first_1;
U_first = utraj_first;

X_second_0T = x_meas_vec_second_0;
X_second_1T = x_meas_vec_second_1;
U_second = utraj_second;

% Plot simulated trajectory (for debugging)
figure;
subplot(1,2,1); hold on; box on; plot(x(1,:),x(2,:),'b'); xlabel('x_1'); ylabel('x_2');
close;

steps_first = length(x_first);
steps_second = length(x_second);

% Construct matrix zonotope for noise (first subsystem)
index = 1;
for i = 1:size(W.generators, 2)
    vec = W.Z(:, i + 1);
    GW1{index} = [vec, zeros(dim_x, steps_first - 1)];
    for j = 1:steps_first - 1
        GW1{j + index} = [GW1{index + j - 1}(:, 2:end) GW1{index + j - 1}(:, 1)];
    end
    index = j + index + 1;
end
Wmatzono1 = matZonotope(zeros(dim_x, steps_first), GW1);

% Construct matrix zonotope for noise (second subsystem)
index = 1;
for i = 1:size(W.generators, 2)
    vec = W.Z(:, i + 1);
    GW2{index} = [vec, zeros(dim_x, steps_second - 1)];
    for j = 1:steps_second - 1
        GW2{j + index} = [GW2{index + j - 1}(:, 2:end) GW2{index + j - 1}(:, 1)];
    end
    index = j + index + 1;
end
Wmatzono2 = matZonotope(zeros(dim_x, steps_second), GW2);

% Compute the system matrices
X1W_cen_1 = X_first_1T - Wmatzono1.center;
X1W_1 = matZonotope(X1W_cen_1, Wmatzono1.generator);

X1W_cen_2 = X_second_1T - Wmatzono2.center;
X1W_2 = matZonotope(X1W_cen_2, Wmatzono2.generator);

% Compute over-approximations of [A, B] matrices
AB1 = X1W_1 * pinv([X_first_0T; U_first]);
AB2 = X1W_2 * pinv([X_second_0T; U_second]);

% Validate that true system matrices are within the approximations
intAB11 = intervalMatrix(AB1);
intAB1 = intAB11.int;
intAB1.sup >= [A1, B1];
intAB1.inf <= [A1, B1];

intAB22 = intervalMatrix(AB2);
intAB2 = intAB22.int;
intAB2.sup >= [A2, B2];
intAB2.inf <= [A2, B2];
end


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



