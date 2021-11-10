% three_pool_cycle_length.m
% Function used to simulate realizations of stochastic three-pool model from
% Barendregt & Thomas, 2021, as given by Eq. (5.1)-(5.3) for a single cycle.

function t = three_pool_cycle_length(gamma,mu,tau_a,Omega)
% Define stoichiometric matrix for the three-pool system:
v_11 = [-1; 1; 0; 0; 0; 0];
v_12 = [-1; 1; 0; 0; 0; 0];
v_13 = [1; -1; 0; 0; 0; 0];
v_21 = [0; 0; -1; 1; 0; 0];
v_22 = [0; 0; -1; 1; 0; 0];
v_23 = [0; 0; 1; -1; 0; 0];
v_31 = [0; 0; 0; 0; -1; 1];
v_32 = [0; 0; 0; 0; -1; 1];
v_33 = [0; 0; 0; 0; 1; -1];
V = [v_11 v_12 v_13 v_21 v_22 v_23 v_31 v_32 v_33];
% Define microscopic rate constants c_j:
c(1) = 1/(tau_a*Omega); c(2) = mu/(tau_a); c(3) = gamma/(tau_a*Omega);
% Initialize simulation:
A = [0; Omega; Omega; 0; Omega; 0];
t = 0;
cycle_check = 0;

while cycle_check ~= 3
    % Calculate reaction hazards:
    a(1) = c(1)*A(1)*A(2);
    a(2) = c(2)*A(1);
    a(3) = c(3)*A(2)*A(4);
    a(4) = c(1)*A(3)*A(4);
    a(5) = c(2)*A(3);
    a(6) = c(3)*A(4)*A(6);
    a(7) = c(1)*A(5)*A(6);
    a(8) = c(2)*A(5);
    a(9) = c(3)*A(6)*A(2);
    % Draw time of next reaction from exponential distribution:
    asum = sum(a);
    t = t+log(1/rand)/asum;
    % Convert reaction hazards into reaction probabilities and draw
    % reaction to occur:
    n = rand; k = 1; Phi = a(1)/asum;
    while Phi < n
        k = k+1;
        Phi = Phi+a(k)/asum;
    end
    % Update state vector given drawn reaction:
    A = A+V(:,k);
    % Check progress along a cycle using threshold planes that trisect the
    % cubic lattice:
    if (A(2)<=A(6)) & (A(2)<A(4)) & (cycle_check == 0)
        cycle_check = 1;
    elseif (A(4)<A(6)) & (A(2)>=A(4)) & (cycle_check == 1)
        cycle_check = 2;
    elseif (A(2)>A(6)) & (A(4)>=A(6)) & (cycle_check == 2)
        cycle_check = 3;
    elseif (A(4)<A(6)) & (A(2)>=A(4)) & (cycle_check == 0)
        cycle_check = -1;
    elseif (A(2)>A(6)) & (A(4)>=A(6)) & (cycle_check == -1)
        cycle_check = 0;
    end
end