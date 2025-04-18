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
