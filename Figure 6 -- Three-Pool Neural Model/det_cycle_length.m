% det_cycle_length.m
% Function used to approximate cycle length of deterministic three-pool
% model given by equation (2.4) of Barendregt & Thomas, 2021.

function T = det_cycle_length(gamma,tau,mu)
% Define initial condition and threshold tolerance:
init = [0.9;0.1;0.1]; tol = 0.9;
for i = 1:length(mu)
    % Simulate trajectory of three-pool model:
    [t,three_pool] = ode45(@(t,three_pool) three_pool_MF(t,three_pool,tau,gamma,mu(i)),[0 100],init);
    t_i = [];
    for j = 1:(length(t)-1)
        % Find times trajectory crosses threshold:
        if three_pool(j,1)<tol && three_pool(j+1,1)>=tol
            t_i = [t_i t(j)];
        end
    end
    % Calculate approximate cycle length:
    T(i) = mean(diff(t_i));
end