% GV_MF.m
% Function used to simulate mean-field of GV model from Barendregt &
% Thomas, 2021, as given in Eq. (2.2).

function dNdt = GV_MF(t,N,alpha,beta,b,d,Omega)
dNdt = [N(1)*(Omega*(b-d)-N(1)-alpha*N(2)-beta*N(3))/Omega;...
    N(2)*(Omega*(b-d)-beta*N(1)-N(2)-alpha*N(3))/Omega;...
    N(3)*(Omega*(b-d)-alpha*N(1)-beta*N(2)-N(3))/Omega];
end