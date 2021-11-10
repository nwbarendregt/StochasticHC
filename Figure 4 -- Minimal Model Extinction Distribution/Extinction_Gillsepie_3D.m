% Extinction_Gillespie_3D.m
% Generate empirical extinction distribution data for three-dimensional 
% minimal model, as shown in Fig. 4C of Barendregt & Thomas, 2021.

clear

% Define model parameters and sample size:
vol = 30; alpha = 0.8; beta = 1.3; r = 1; Nint = [vol; vol; vol]/3;
nsamp = 1e6; 
% Pre-allocate empirical distribution storage:
pi_x = zeros(2*vol+1); pi_y = pi_x; pi_z = pi_x;
% Set simulation to calculate 0 extinction -> 1 extinction distribution:
Next = 0;

for i = 1:nsamp
    % Simulate realization of minimal model:
    [~,N] = minimal_gillespie_extinction(alpha,beta,r,vol,Nint,Next);
    % Store simulation in appropriate distribution:
    if N(1) == 0
        pi_x(N(2)+1,N(3)+1) = pi_x(N(2)+1,N(3)+1)+1/nsamp;
    elseif N(2) == 0
        pi_y(N(1)+1,N(3)+1) = pi_y(N(1)+1,N(3)+1)+1/nsamp;
    else
        pi_z(N(1)+1,N(2)+1) = pi_z(N(1)+1,N(2)+1)+1/nsamp;
    end
    disp(i)
end

save('Extinction_Gillespie_3D_Data')