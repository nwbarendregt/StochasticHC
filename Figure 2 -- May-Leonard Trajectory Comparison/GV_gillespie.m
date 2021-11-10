% GV_gillespie.m
% Function use to simulate realizations of stochastic GV model from
% Barendregt & Thomas, 2021, as given by Eq. (3.1)-(3.5).

function [t,N] = GV_gillespie(alpha,beta,b,d,vol,Nint,t_f)
% Define stoichiometric matrix for the GV system:
v_1b = [1; 0; 0];
v_1d = [-1; 0; 0];
v_11 = [-1; 0; 0];
v_12 = [-1; 0; 0];
v_13 = [-1; 0; 0];
v_2b = [0; 1; 0];
v_2d = [0; -1; 0];
v_21 = [0; -1; 0];
v_22 = [0; -1; 0];
v_23 = [0; -1; 0];
v_3b = [0; 0; 1];
v_3d = [0; 0; -1];
v_31 = [0; 0; -1];
v_32 = [0; 0; -1];
v_33 = [0; 0; -1];
V = [v_1b v_1d v_11 v_12 v_13 v_2b v_2d v_21 v_22 v_23 v_3b v_3d v_31 v_32 v_33];

% Define microscopic rate constants c_j:
c(1) = b; c(2) = d; c(3) = 2/vol; c(4) = alpha/vol; c(5) = beta/vol;
c(6) = b; c(7) = d; c(8) = beta/vol; c(9) = 2/vol; c(10) = alpha/vol;
c(11) = b; c(12) = d; c(13) = alpha/vol; c(14) = beta/vol; c(15) = 2/vol;

% Initialize simulation:
t = 0;
N = Nint;
i = 1;

while t(i) < t_f
    i = i+1;
    % Calculate reaction hazards:
    a(1) = c(1)*N(1,i-1);
    a(2) = c(2)*N(1,i-1);
    a(3) = c(3)*N(1,i-1)*(N(1,i-1)-1)/2;
    a(4) = c(4)*N(1,i-1)*N(2,i-1);
    a(5) = c(5)*N(1,i-1)*N(3,i-1);
    a(6) = c(6)*N(2,i-1);
    a(7) = c(7)*N(2,i-1);
    a(8) = c(8)*N(2,i-1)*N(1,i-1);
    a(9) = c(9)*N(2,i-1)*(N(2,i-1)-1)/2;
    a(10) = c(10)*N(2,i-1)*N(3,i-1);
    a(11) = c(11)*N(3,i-1);
    a(12) = c(12)*N(3,i-1);
    a(13) = c(13)*N(3,i-1)*N(1,i-1);
    a(14) = c(14)*N(3,i-1)*N(2,i-1);
    a(15) = c(15)*N(3,i-1)*(N(3,i-1)-1)/2;
    % Draw time of next reaction from exponential distribution:
    asum = sum(a);
    t(i) = t(i-1)+log(1/rand)/asum;
    % Convert reaction hazards into reaction probabilities and draw
    % reaction to occur:
    n = rand; k = 1; Phi = a(1)/asum;
    while Phi < n
        k = k+1;
        Phi = Phi+a(k)/asum;
    end
    % Update state vector given drawn reaction:
    N(:,i) = N(:,i-1)+V(:,k);
end