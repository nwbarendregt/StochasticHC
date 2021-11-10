% ML_MF.m
% Function used to simulate mean-field of May-Leonard system from 
% Barendregt & Thomas, 2021, as given in Eq. (1.1).

function dNdt = ML_MF(t,N,alpha,beta)
dNdt = [N(1)*(1-N(1)-alpha*N(2)-beta*N(3));...
    N(2)*(1-beta*N(1)-N(2)-alpha*N(3));...
    N(3)*(1-alpha*N(1)-beta*N(2)-N(3))];
end