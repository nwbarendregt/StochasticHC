% Figure_5_Generate.m
% Script used to generate Figure 5 from Barendregt & Thomas, 2021. 

clear
% Load data for Fig. 5A:
load('Extinction_Times_Discrete_3D_Data.mat')
% Generate Fig. 5A:
figure
mesh(X,Y,Z,C,'edgecolor','interp','facecolor','interp')
xlim([1,30]); ylim([1,30]); zlim([1,30])
colorbar
view(115,30)
% Load data for Fig. 5B:
load('Extinction_Times_Gillespie_3D_Data.mat')
% Generate Fig. 5B:
figure
mesh(X,Y,Z,mean(C,3),'edgecolor','interp','facecolor','interp')
xlim([1,30]); ylim([1,30]); zlim([1,30])
colorbar
view(115,30)
% Load data for Fig. 5C:
load('Extinction_Times_Discrete_3D_Data.mat','C'); C_D = C;
load('Extinction_Times_Gillespie_3D_Data.mat','C'); C_G = mean(C,3);
% Generate Fig. 5C:
figure
scatter(C_D(:),C_G(:),'filled')
hold on
plot([0 8],[0 8],'--')
% Load data for Fig. 5D:
load('Extinction_Times_Gillespie_vol_Data.mat')
% Generate Fig. 5D:
figure
loglog(vol,T_avg)
xlim([1e1 1e5])