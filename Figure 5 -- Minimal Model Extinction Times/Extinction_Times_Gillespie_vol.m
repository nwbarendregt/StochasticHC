% Extinction_Times_Gillespie_vol.m
% Script used to generate data for Fig. 5D in Barendregt & Thomas, 2021. 

clear
% Define simulation parameters:
alpha = 0.8; beta = 1.3; r = 1;
nsamp = 1e4; Next = 0;
vol = unique(12*round(logspace(0,4))); s = [0.25 0.5 0.75 1];
% Pre-allocate storage of extinction time matrix:
T = NaN(length(s),length(vol),nsamp);

for k = 1:nsamp
    T_k = NaN(length(s),length(vol));
    for i = 1:length(s)
        for j = 1:length(vol)
            % Compute initial condition along segment described in text:
            x_init = vol(j)*(1-s(i))+vol(j)*s(i)/3;
            y_init = vol(j)*s(i)/3;
            z_init = y_init;
            % Simulate realization of minimal model:
            T_k(i,j) = minimal_gillespie_extinction(alpha,beta,r,vol(j),[x_init; y_init; z_init],Next);
        end
    end
    T(:,:,k) = T_k;
end
% Compute average extinction times:
T_avg = mean(T,3);
save('Extinction_Times_Gillespie_vol_Data');