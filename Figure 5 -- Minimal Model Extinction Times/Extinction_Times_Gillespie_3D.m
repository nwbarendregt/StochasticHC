% Extinction_Times_Gillespie_3D.m
% Script used to generate data for Fig. 5B,C in Barendregt & Thomas, 2021.

clear
% Define simulation parameters:
vol = 30; alpha = 0.8; beta = 1.3; r = 1;
nsamp = 1e4; Next = 0;
% Find initial state vectors on plane Pi as described in text:
[x,y,z] = meshgrid(0:(2*vol),0:(2*vol),0:(2*vol));
x = x(:); y = y(:); z = z(:);
X_int = []; Y_int = []; Z_int = [];
for i = 1:length(x)
    if x(i)+y(i)+z(i) == vol
        X_int = [X_int x(i)];
        Y_int = [Y_int y(i)];
        Z_int = [Z_int z(i)];
    end
end
% Pre-allocate storage of extinction time matrix:
T = NaN(length(X_int),nsamp);
for i = 1:length(X_int)
    for j = 1:nsamp
        % Simulate minimal model:
        T(i,j) = minimal_gillespie_extinction(alpha,beta,r,vol,[X_int(i); Y_int(i); Z_int(i)],Next);
    end
end 
% Project slice of empirical extinction times onto the plane Pi:
[X,Y] = meshgrid(0:(2*vol),0:(2*vol)); Z = NaN(2*vol+1,2*vol+1); C = NaN([size(Z),nsamp]); k = 1;
while k <= (length(X_int))
    for ii = 1:vol
        for jj = 1:vol
            if (X(ii,jj) == X_int(k)) && (Y(ii,jj) == Y_int(k))
                Z(ii,jj) = Z_int(k); C(ii,jj,:) = T(k,:);
            end
        end
    end
    k = k+1;
end 

save('Extinction_Times_Gillespie_3D_Data')