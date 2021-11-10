% Fig_1_Generate.m
% Script used to generate Figure 1 from Barendregt & Thomas, 2021.

clear
% Define model parameters and initial condition:
alpha = 0.8; beta = 1.3; n_init = [1;0.8;0.2];
% Simulate mean-field equations for May-Leonard system:
[t,ML] = ode45(@(t,ML) ML_MF(t,ML,alpha,beta),[0 200],n_init);
% Plot Fig. 1A:
figure
plot(t,ML)

% Define model parameters and initial condition:
tau = 1; gamma = 2.4; mu = 1e-5; a_init = n_init;
% Simulate mean-field equations for three-pool model:
[t,three_pool] = ode45(@(t,CPG) CPG_MF(t,CPG,tau,gamma,mu),[0 200],a_init);
% Plot Fig. 1B:
figure
plot(t,three_pool)

% Plot Fig. 1C:
figure
plot3(ML(:,1),ML(:,2),ML(:,3))
view(115,30)
grid on
% Plot Fig. 1D:
figure
plot3(three_pool(:,1),three_pool(:,2),three_pool(:,3))
view(115,30)
grid on