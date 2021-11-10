% Figure_4_Generate.m
% Script used to generate Figure 4 from Barendregt & Thomas, 2021.

% Load data for Fig. 4A:
load('Extinction_Distribution_Discrete_3D_Data.mat')
% Plot Fig. 4A:
figure
[X,Y,Z] = meshgrid(0:(2*vol-1),0:(2*vol-1),0:(2*vol-1));
Pi = NaN(2*vol,2*vol,2*vol);
Pi(:,1,:) = squeeze(P_abs(1,2:end,2:end));
xslice = 0;
h = slice(X,Y,Z,Pi,xslice,[],[]);
set(h,'edgecolor','none')
hold on
Pi = NaN(2*vol,2*vol,2*vol);
Pi(1,:,:) = squeeze(P_abs(2:end,1,2:end));
yslice = 0;
h = slice(X,Y,Z,Pi,[],yslice,[]);
set(h,'edgecolor','none')
Pi = NaN(2*vol,2*vol,2*vol);
Pi(:,:,1) = P_abs(2:end,2:end,1)';
zslice = 0;
h = slice(X,Y,Z,Pi,[],[],zslice);
set(h,'edgecolor','none')
axis([0 vol 0 vol 0 vol])
colorbar
xlabel('N_1'); ylabel('N_2'); zlabel('N_3');
view(115,30)

% Load data for Fig. 4B:
load('Extinction_Distribution_Discrete_2D_Data.mat')
% Plot Fig. 4B:
figure
plot(P_abs(:,1))
hold on
plot(P_abs(1,:))

% Load data for Fig. 4C:
load('Extinction_Gillespie_3D_Data.mat')
% Plot Fig. 4C:
figure
[X,Y,Z] = meshgrid(0:(2*vol-1),0:(2*vol-1),0:(2*vol-1));
Pi = NaN(2*vol,2*vol,2*vol);
Pi(:,1,:) = pi_x(2:end,2:end)';
xslice = 0;
h = slice(X,Y,Z,Pi,xslice,[],[]);
set(h,'edgecolor','none')
hold on
Pi = NaN(2*vol,2*vol,2*vol);
Pi(1,:,:) = pi_y(2:end,2:end)';
yslice = 0;
h = slice(X,Y,Z,Pi,[],yslice,[]);
set(h,'edgecolor','none')
Pi = NaN(2*vol,2*vol,2*vol);
Pi(:,:,1) = pi_z(2:end,2:end);
zslice = 0;
h = slice(X,Y,Z,Pi,[],[],zslice);
set(h,'edgecolor','none')
axis([0 vol 0 vol 0 vol])
colorbar
xlabel('N_1'); ylabel('N_2'); zlabel('N_3');
view(115,30)

% Load data for Fig. 4D:
load('Extinction_Gillespie_2D_Data.mat')
% Plot Fig. 4D:
P = histcounts2(pi,xy,-0.5:(max(pi)+0.5),0.5:2.5,'normalization','probability');
figure
plot(P(:,1))
hold on
plot(P(:,2))
xlim([0 2*vol])