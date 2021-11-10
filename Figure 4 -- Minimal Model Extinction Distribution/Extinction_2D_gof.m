% Extinction_2D_gof.m
% Script used to calculate the goodness-of-fit between analytic
% distribution in Fig. 4B and empirical distribution in Fig. 4D of 
% Barendregt & Thomas, 2021.

clear
% Load empirical distribution:
load('Extinction_Gillespie_2D_Data.mat');
% Calculate observed frequency of each state's extinction:
P = histcounts2(pi,xy,-0.5:(2*vol+0.5),0.5:2.5);
% Load analytic distribution:
load('Extinction_Distribution_Discrete_2D_Data.mat','P_abs');
% Calculate expected frequency of each state's extinction:
pi_x = nsamp*squeeze(P_abs(:,1)); pi_y = nsamp*squeeze(P_abs(1,:));
% Calculate p-values from chi-squared goodness-of-fit test:
[~,~,stx] = chi2gof(1:(2*vol+1),'ctrs',1:(2*vol+1),'frequency',P(:,1)',...
        'expected',pi_x(:)');
[~,~,sty] = chi2gof(1:(2*vol+1),'ctrs',1:(2*vol+1),'frequency',P(:,2)',...
        'expected',pi_y(:)');

obsCounts = [stx.O sty.O]; expCounts = [stx.E sty.E]; bins = 1:length(obsCounts);
[~,p,~] = chi2gof(bins,'ctrs',bins,'frequency',obsCounts,'expected',expCounts);
disp(p)