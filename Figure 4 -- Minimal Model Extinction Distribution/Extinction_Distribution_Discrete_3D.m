% Extinction_Distribution_Discrete_3D.m
% Script used to construct and solve the first-passage location problem in
% Barendregt & Thomas, 2021, as given by Eq. (4.1).
% Infinitesimal generator matrix L for the full three-dimensional system 
% constructed using the approach described in Appendix C.

clear
% Define model parameters and initial state at deterministic fixed-point:
vol = 30; alpha = 0.8; beta = 1.3; r = 1;
init_state = vol/3*[1 1 1];
% Calculate microscopic rate constants c_j:
c(1) = r; c(2) = 2/vol; c(3) = alpha/vol; c(4) = beta/vol;
c(5) = r; c(6) = beta/vol; c(7) = 2/vol; c(8) = alpha/vol;
c(9) = r; c(10) = alpha/vol; c(11) = beta/vol; c(12) = 2/vol;
% Define simulation domain [0, 2*Omega]^3:
N1 = 0:(2*vol); N2 = N1; N3 = N2;
[N1_m,N2_m,N3_m] = meshgrid(N1,N2,N3);
% Generate state vector and find index of initial state:
states = [N1_m(:) N2_m(:) N3_m(:)];
init_i = 1+(1+2*vol)^2*init_state(3)+(1+2*vol)*init_state(1)+init_state(2);

% Pre-allocate storage of infinitesimal generator matrix A:
A = sparse(length(states),length(states));
% Construct A row-by-row:
for i = 1:length(states)
    A_i = zeros(1,length(states));
    if sum(states(i,:)==0)~=0 % State i is on absorbing boundary
        A_i(i) = 1;
    elseif sum(states(i,:)==2*vol)==3 % State i is 2*Omega*[1,1,1]
        A_i(i-(1+2*vol)^2) = c(10)*states(i,3)*states(i,1)+...
            c(11)*states(i,3)*states(i,2)+c(12)*states(i,3)*(states(i,3)-1)/2;
        A_i(i-(1+2*vol)) = c(2)*states(i,1)*(states(i,1)-1)/2+...
            c(3)*states(i,1)*states(i,2)+c(4)*states(i,1)*states(i,3);
        A_i(i-1) = c(6)*states(i,2)*states(i,1)+...
            c(7)*states(i,2)*(states(i,2)-1)/2+c(8)*states(i,2)*states(i,3);
        A_i(i) = -sum(A_i);
    elseif sum(states(i,:)==2*vol)==2
        if states(i,1)~=2*vol % State i is on intersection of N_2 = 2*Omega and N_3 = 2*Omega
            A_i(i-(1+2*vol)^2) = c(10)*states(i,3)*states(i,1)+...
                c(11)*states(i,3)*states(i,2)+c(12)*states(i,3)*(states(i,3)-1)/2;
            A_i(i-(1+2*vol)) = c(2)*states(i,1)*(states(i,1)-1)/2+...
                c(3)*states(i,1)*states(i,2)+c(4)*states(i,1)*states(i,3);
            A_i(i-1) = c(6)*states(i,2)*states(i,1)+...
                c(7)*states(i,2)*(states(i,2)-1)/2+c(8)*states(i,2)*states(i,3);
            A_i(i+(1+2*vol)) = c(1)*states(i,1);
            A_i(i) = -sum(A_i);
        elseif states(i,2)~=2*vol % State i is on intersection of N_1 = 2*Omega and N_3 = 2*Omega
            A_i(i-(1+2*vol)^2) = c(10)*states(i,3)*states(i,1)+...
                c(11)*states(i,3)*states(i,2)+c(12)*states(i,3)*(states(i,3)-1)/2;
            A_i(i-(1+2*vol)) = c(2)*states(i,1)*(states(i,1)-1)/2+...
                c(3)*states(i,1)*states(i,2)+c(4)*states(i,1)*states(i,3);
            A_i(i-1) = c(6)*states(i,2)*states(i,1)+...
                c(7)*states(i,2)*(states(i,2)-1)/2+c(8)*states(i,2)*states(i,3);
            A_i(i+1) = c(5)*states(i,2);
            A_i(i) = -sum(A_i);
        else % State i is on intersection of N_1 = 2*Omega and N_2 = 2*Omega
            A_i(i-(1+2*vol)^2) = c(10)*states(i,3)*states(i,1)+...
                c(11)*states(i,3)*states(i,2)+c(12)*states(i,3)*(states(i,3)-1)/2;
            A_i(i-(1+2*vol)) = c(2)*states(i,1)*(states(i,1)-1)/2+...
                c(3)*states(i,1)*states(i,2)+c(4)*states(i,1)*states(i,3);
            A_i(i-1) = c(6)*states(i,2)*states(i,1)+...
                c(7)*states(i,2)*(states(i,2)-1)/2+c(8)*states(i,2)*states(i,3);
            A_i(i+(1+2*vol)^2) = c(9)*states(i,3);
            A_i(i) = -sum(A_i);
        end
    elseif sum(states(i,:)==2*vol)==1
        if states(i,1)==2*vol % State i is on reflecting boundary N_1 = 2*Omega
            A_i(i-(1+2*vol)^2) = c(10)*states(i,3)*states(i,1)+...
                c(11)*states(i,3)*states(i,2)+c(12)*states(i,3)*(states(i,3)-1)/2;
            A_i(i-(1+2*vol)) = c(2)*states(i,1)*(states(i,1)-1)/2+...
                c(3)*states(i,1)*states(i,2)+c(4)*states(i,1)*states(i,3);
            A_i(i-1) = c(6)*states(i,2)*states(i,1)+...
                c(7)*states(i,2)*(states(i,2)-1)/2+c(8)*states(i,2)*states(i,3);
            A_i(i+1) = c(5)*states(i,2);
            A_i(i+(1+2*vol)^2) = c(9)*states(i,3);
            A_i(i) = -sum(A_i);
        elseif states(i,2)==2*vol % State i is on reflecting boundary N_2 = 2*Omega
            A_i(i-(1+2*vol)^2) = c(10)*states(i,3)*states(i,1)+...
                c(11)*states(i,3)*states(i,2)+c(12)*states(i,3)*(states(i,3)-1)/2;
            A_i(i-(1+2*vol)) = c(2)*states(i,1)*(states(i,1)-1)/2+...
                c(3)*states(i,1)*states(i,2)+c(4)*states(i,1)*states(i,3);
            A_i(i-1) = c(6)*states(i,2)*states(i,1)+...
                c(7)*states(i,2)*(states(i,2)-1)/2+c(8)*states(i,2)*states(i,3);
            A_i(i+(1+2*vol)) = c(1)*states(i,1);
            A_i(i+(1+2*vol)^2) = c(9)*states(i,3);
            A_i(i) = -sum(A_i);
        else % State i is on reflecting boundary N_3 = 2*Omega
            A_i(i-(1+2*vol)^2) = c(10)*states(i,3)*states(i,1)+...
                c(11)*states(i,3)*states(i,2)+c(12)*states(i,3)*(states(i,3)-1)/2;
            A_i(i-(1+2*vol)) = c(2)*states(i,1)*(states(i,1)-1)/2+...
                c(3)*states(i,1)*states(i,2)+c(4)*states(i,1)*states(i,3);
            A_i(i-1) = c(6)*states(i,2)*states(i,1)+...
                c(7)*states(i,2)*(states(i,2)-1)/2+c(8)*states(i,2)*states(i,3);
            A_i(i+1) = c(5)*states(i,2);
            A_i(i+(1+2*vol)) = c(1)*states(i,1);
            A_i(i) = -sum(A_i);
        end
    else % State i is on the interior of the domain
        A_i(i-(1+2*vol)^2) = c(10)*states(i,3)*states(i,1)+...
            c(11)*states(i,3)*states(i,2)+c(12)*states(i,3)*(states(i,3)-1)/2;
        A_i(i-(1+2*vol)) = c(2)*states(i,1)*(states(i,1)-1)/2+...
            c(3)*states(i,1)*states(i,2)+c(4)*states(i,1)*states(i,3);
        A_i(i-1) = c(6)*states(i,2)*states(i,1)+...
            c(7)*states(i,2)*(states(i,2)-1)/2+c(8)*states(i,2)*states(i,3);
        A_i(i+1) = c(5)*states(i,2);
        A_i(i+(1+2*vol)) = c(1)*states(i,1);
        A_i(i+(1+2*vol)^2) = c(9)*states(i,3);
        A_i(i) = -sum(A_i);
    end
    A(i,:) = A_i;
end

% Find LU decomposition of A to solve first-passage location problem for
% many right-hand sides:
[L,U,P] = lu(A);
% Pre-allocate storage of extinction distribution:
P_abs = NaN(length(states),1);
% Solve first-passage location problem for each absorbing state:
for i = 1:length(states)
    if sum(states(i,:)==0)~=0
        % Construct right-hand side e_i:
        b = sparse(length(states),1); b(i) = 1;
        % Solve first-passage location problem:
        q = L\(P*b); p = U\q; P_abs(i) = p(init_i);
    end
end
clear L U
% Reshape distribution for plotting:
P_abs = reshape(P_abs,[2*vol+1,2*vol+1,2*vol+1]);
% save('Extinction_Distribution_Discrete_3D_Data');