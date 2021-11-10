% Three_Pool_Cycle_Length_Gillespie.m
% Script used to calculate empirical cycle length statistics for three-pool
% model in Barendregt & Thomas, 2021.

clear
% Define simulation parameters:
gamma = 2.4; tau_a = 0.05;
mu = logspace(-8,-3,100); Omega = unique(floor(logspace(0,5,100)));
n_trial = 1e4;
% Pre-allocate storage of cycle length matrix:
t_data = NaN(length(Omega),length(mu),n_trial);
for i = 1:n_trial
    T = NaN(length(Omega),length(mu));
    for j = 1:length(Omega)
        for k = 1:length(mu)
            % Simulate realization of three-pool model:
            T(j,k) = three_pool_cycle_length(gamma,mu(k),tau_a,Omega(j));
        end
    end
    t_data(:,:,i) = T;
end
% Calculate cycle length statistics and delete bulk data for storage:
T_mean = mean(t_data,3); T_var = var(t_data,[],3); T_cv = sqrt(T_var)./T_mean;
clear t_data

save('Three_Pool_Cycle_Length_Data');