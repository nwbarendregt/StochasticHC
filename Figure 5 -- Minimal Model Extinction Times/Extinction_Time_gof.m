% Extinction_Time_gof.m
% Script used to calculate goodness-of-fit between mean extinction times in
% Fig. 5A,B in Barendregt & Thomas, 2021.

clear
% Load mean data:
load('Extinction_Times_Discrete_3D_Data.mat','C'); m = C;
load('Extinction_Times_Gillespie_3D_Data.mat','C'); x = C;
% Pre-allocate p-value storage:
[I,J] = size(m); p = NaN(I,J);

for i = 1:I
    for j = 1:J
        if ~isnan(m(i,j))
            % Compute p-value using t-test:
            [~,p(i,j)] = ttest(x(i,j,:),m(i,j));
        end
    end
end
% Compute average p-value over the plane Pi:
p = mean(p,'all','omitnan');