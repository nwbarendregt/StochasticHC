% Extinction_Times_Discrete_3D.m
% Script used to construct and solve the first-passage time problem in
% Barendregt & Thomas, 2021, as given by Eq. (4.2).
% Infinitesimal generator matrix L for the reduced two-dimensional system
% constructed using the approach described in Appendix C.

clear
% Load infinitesimal generator matrix L:
load('Extinction_Distribution_Discrete_3D_Data.mat');
% Define simulation domain:
[x,y,z] = meshgrid(0:(2*vol),0:(2*vol),0:(2*vol));
x = x(:); y = y(:); z = z(:);
% Find initial state vectors on plane Pi as described in text:
X_int = []; Y_int = []; Z_int = [];
for i = 1:length(x)
    if x(i)+y(i)+z(i) == vol
        X_int = [X_int x(i)];
        Y_int = [Y_int y(i)];
        Z_int = [Z_int z(i)];
    end
end
% Construct right-hand side b:
b = zeros(length(states),1);
for i = 1:length(states)
    if sum(states(i,:)==0)==0
        b(i) = -1;
    end
end
% Solve first-passage time problem:
T = A\b;
% Take slice of solution to first-passage time problem along the plane Pi:
[X,Y] = meshgrid(0:(2*vol),0:(2*vol)); Z = NaN(2*vol+1,2*vol+1); C = Z; k = 1;
while k <= (length(X_int))
    for ii = 1:length(X)
        for jj = 1:length(Y)
            if (X(ii,jj) == X_int(k)) && (Y(ii,jj) == Y_int(k))
                Z(ii,jj) = Z_int(k);
                i = find(sum(states == [X_int(k) Y_int(k) Z_int(k)],2)==3);
                C(ii,jj) = T(i);
            end
        end
    end
    k = k+1;
end

save('Extinction_Times_Discrete_3D_Data')