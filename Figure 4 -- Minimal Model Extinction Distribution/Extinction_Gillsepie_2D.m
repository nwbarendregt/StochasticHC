% Extinction_Gillespie_3D.m
% Generate empirical extinction distribution data for reduced
% two-dimensional minimal model with N_3 going extinct first, as shown in
% Fig. 4D of Barendregt & Thomas, 2021.

clear
% Load empirical three-dimensional distribution to calculate likelihood of
% each initial state conditioned on N_3 going extinct first:
load('Extinction_Gillespie_3D_Data.mat'); pi_z = pi_z(:)/sum(pi_z,'all');
% Pre-allocate empirical distribution storage:
pi = NaN(1,nsamp); xy = NaN(1,nsamp);
% Set simulation to calculate 1 extinction -> 2 extinctions distribution:
Next = 1;

for i = 1:nsamp
    % Sample initial state to seed simulation from:
    n = rand; k = 1; Phi = pi_z(k);
    while Phi < n
        k = k+1;
        Phi = Phi+pi_z(k);
    end
    [int1,int2] = ind2sub([2*vol+1,2*vol+2],k);
    % Simulate realization of minimal model:
    [~,N] = minimal_gillespie_extinction(alpha,beta,r,vol,[int1-1;int2-1;0],Next);
    % Store simulation in appropriate distribution:
    if N(1) == 0
        xy(i) = 1; pi(i) = N(2);
    else
        xy(i) = 2; pi(i) = N(1);
    end
    disp(i)
end

save('Extinction_Gillespie_2D_Data')