% Fig_2_Generate.m
% Script used to generate Figure 2 from Barendregt & Thomas, 2021.
% Script produces panels A, C, and E of figure; to obtain panels B, D, and
% F, uncomment line 16, 23, and 30, respectively.

clear 
% Define model parameters:
alpha = 0.8; beta = 1.3; b = 3; d = 2; Omega = 30;
% Set simulation time and initial condition:
t_f = 200; Nint = [1/3;1/3;1/3+0.1]*Omega;
% Simulate mean-field equations for GV model:
[t,N] = ode45(@(t,X) GV_MF(t,X,alpha,beta,b,d,Omega),[0 t_f],Nint);
% Plot Fig. 2A:
figure
plot(t,N)
% xlim([0 20]) % Uncomment line to generate Fig. 2B.

% Simulate realization of stochastic GV system:
[T,N] = GV_Gillespie(alpha,beta,b,d,Omega,Nint,t_f);
% Plot Fig. 2C:
figure
plot(T,N)
% xlim([0 20]) % Uncomment line to generate Fig. 2D.

% Simulate realization of stochastic minimal system:
[T,N] = GV_Gillespie(alpha,beta,b-d,0,Omega,Nint,t_f);
% Plot Fig. 2E:
figure
plot(T,N)
% xlim([0 20]) % Uncomment line to generate Fig. 2F.