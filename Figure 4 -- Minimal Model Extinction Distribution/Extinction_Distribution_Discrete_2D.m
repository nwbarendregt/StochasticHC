% Extinction_Distribution_Discrete_2D.m
% Script used to construct and solve the first-passage location problem in
% Barendregt & Thomas, 2021, as given by Eq. (4.1).
% Infinitesimal generator matrix L for the reduced two-dimensional system
% constructed using the approach described in Appendix C.

clear
% Define model parameters:
vol = 30; alpha = 0.8; beta = 1.3; r = 1;
% Calculate microscopic rate constants c_j:
c(1) = r; c(2) = 2/vol; c(3) = alpha/vol;
c(5) = r; c(6) = beta/vol; c(7) = 2/vol;
% Define simulation domain [0, 2*Omega]^2:
N1 = 0:(2*vol); N2 = N1;
[N1_m,N2_m] = meshgrid(N1,N2);
% Generate state vector:
states = [N1_m(:) N2_m(:)];
% Load full three-dimensional distribution to calculate likelihood of each
% initial state conditioned on N_3 going extinct first:
load('Extinction_Distribution_Discrete_3D_Data.mat','P_abs')
P_z = P_abs(:,:,1); N = sum(P_z,'all'); P_z = P_z(:);

% Pre-allocate storage of infinitesimal generator matrix A:
A = sparse(length(states),length(states));
% Construct A row-by-row:
for i = 1:length(states)
    A_i = zeros(1,length(states));
    if sum(states(i,:)==0)~=0 % State i is on absorbing boundary
        A_i(i) = 1;
    elseif sum(states(i,:)==2*vol)==2 % State i is 2*Omega*[1,1]
        A_i(i-(1+2*vol)) = c(2)*states(i,1)*(states(i,1)-1)/2+...
            c(3)*states(i,1)*states(i,2);
        A_i(i-1) = c(6)*states(i,2)*states(i,1)+...
            c(7)*states(i,2)*(states(i,2)-1)/2;
        A_i = A_i/sum(A_i);
        A_i(i) = -1;
    elseif sum(states(i,:)==2*vol)==1
        if states(i,1)==2*vol % State i is on reflecting boundary N_1 = 2*Omega
            A_i(i-(1+2*vol)) = c(2)*states(i,1)*(states(i,1)-1)/2+...
                c(3)*states(i,1)*states(i,2);
            A_i(i-1) = c(6)*states(i,2)*states(i,1)+...
                c(7)*states(i,2)*(states(i,2)-1)/2;
            A_i(i+1) = c(5)*states(i,2);
            A_i = A_i/sum(A_i);
            A_i(i) = -1;
        else % State i is on reflecting boundary N_1 = 2*Omega
            A_i(i-(1+2*vol)) = c(2)*states(i,1)*(states(i,1)-1)/2+...
                c(3)*states(i,1)*states(i,2);
            A_i(i-1) = c(6)*states(i,2)*states(i,1)+...
                c(7)*states(i,2)*(states(i,2)-1)/2;
            A_i(i+(1+2*vol)) = c(1)*states(i,1);
            A_i = A_i/sum(A_i);
            A_i(i) = -1;
        end
    else % State i is on the interior of the domain
        A_i(i-(1+2*vol)) = c(2)*states(i,1)*(states(i,1)-1)/2+...
            c(3)*states(i,1)*states(i,2);
        A_i(i-1) = c(6)*states(i,2)*states(i,1)+...
            c(7)*states(i,2)*(states(i,2)-1)/2;
        A_i(i+1) = c(5)*states(i,2);
        A_i(i+(1+2*vol)) = c(1)*states(i,1);
        A_i = A_i/sum(A_i);
        A_i(i) = -1;
    end
    A(i,:) = A_i;
end

% Find LU decomposition of A to solve first-passage location problem for
% many right-hand sides:
[L,U,P] = lu(A);
% Pre-allocate storage of extinction distribution:
P_abs = zeros(length(states),1);
% Solve first-passage location problem for each absorbing state:
for i = 1:length(states)
    if sum(states(i,:)==0)~=0
        % Construct right-hand side e_i:
        b = sparse(length(states),1); b(i) = 1;
        % Solve first-passage location problem:
        q = L\(P*b); p = U\q;
        % Marginalize distribution over likelihood of initial state,
        % conditioned on N_3 going extinct first:
        for j = 1:length(states)
            if P_z(j) ~= 0
                P_abs(i) = P_abs(i)+p(j)*P_z(j);
            end
        end
    end
end
clear L U
% Reshape distribution for plotting and re-normalize:
P_abs = reshape(P_abs,[2*vol+1,2*vol+1]); P_abs = P_abs/N;
save('Extinction_Distribution_Discrete_2D_Data');