% minimal_gillespie_extinction.m
% Function used to simulate realizations of stochastic minimal model from
% Barendregt & Thomas, 2021, as given by Eq. (3.6)-(3.9) until extinction.

function [t,N] = minimal_gillespie_extinction(alpha,beta,r,vol,Nint,Next)
% Define stoichiometric matrix for the minimal system:
v_1b = [1; 0; 0];
v_11 = [-1; 0; 0];
v_12 = [-1; 0; 0];
v_13 = [-1; 0; 0];
v_2b = [0; 1; 0];
v_21 = [0; -1; 0];
v_22 = [0; -1; 0];
v_23 = [0; -1; 0];
v_3b = [0; 0; 1];
v_31 = [0; 0; -1];
v_32 = [0; 0; -1];
v_33 = [0; 0; -1];
V = [v_1b v_11 v_12 v_13 v_2b v_21 v_22 v_23 v_3b v_31 v_32 v_33];

% Define microscopic rate constants c_j:
c(1) = r; c(2) = 2/vol; c(3) = alpha/vol; c(4) = beta/vol;
c(5) = r; c(6) = beta/vol; c(7) = 2/vol; c(8) = alpha/vol;
c(9) = r; c(10) = alpha/vol; c(11) = beta/vol; c(12) = 2/vol;

% Initialize simulation:
t = 0;
N = Nint;

while sum(N==0) <= Next
    % Calculate reaction hazards:
    a(1) = c(1)*N(1);
    a(2) = c(2)*N(1)*(N(1)-1)/2;
    a(3) = c(3)*N(1)*N(2);
    a(4) = c(4)*N(1)*N(3);
    a(5) = c(5)*N(2);
    a(6) = c(6)*N(2)*N(1);
    a(7) = c(7)*N(2)*(N(2)-1)/2;
    a(8) = c(8)*N(2)*N(3);
    a(9) = c(9)*N(3);
    a(10) = c(10)*N(3)*N(1);
    a(11) = c(11)*N(3)*N(2);
    a(12) = c(12)*N(3)*(N(3)-1)/2;
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
    N = N+V(:,k);
end