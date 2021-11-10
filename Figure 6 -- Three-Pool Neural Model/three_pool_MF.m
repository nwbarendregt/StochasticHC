% three_pool_MF.m
% Function used to simulate mean-field of three-pool model from Barendregt 
% & Thomas, 2021, as given in Eq. (2.4).

function dAdt = three_pool_MF(t,A,tau,gamma,mu)
dAdt = [(A(1)*(1-A(1)-gamma*A(2))+mu*(1-A(1)))/tau;...
    (A(2)*(1-A(2)-gamma*A(3))+mu*(1-A(2)))/tau;...
    (A(3)*(1-A(3)-gamma*A(1))+mu*(1-A(3)))/tau];
end