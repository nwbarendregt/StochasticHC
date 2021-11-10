% Figure_6_Generate.m
% Script used to generate Figure 6 from Barendregt & Thomas, 2021.

clear
% Load data:
load('Three_Pool_Cycle_Length_Data.mat')
% Generate Fig. 6A:
figure
surf(log10(Omega),log10(mu),log10(T_mean'));
xlim(log10([min(Omega) max(Omega)]))
view(100,20)
% Generate Fig. 6B:
figure
surf(log10(Omega),log10(mu),log10(T_var'));
xlim(log10([min(Omega) max(Omega)]))
view(100,20)
% Generate Fig. 6C:
figure
loglog(Omega,T_mean(:,1),'linewidth',15)
hold on
loglog(Omega,T_mean(:,51),'linewidth',15)
loglog(Omega,T_mean(:,end),'linewidth',15)
loglog(Omega,3*tau_a./(Omega*mu(1)))
loglog(Omega,3*tau_a./(Omega*mu(51)))
loglog(Omega,3*tau_a./(Omega*mu(end)))
loglog(Omega,det_cycle_length(gamma,tau_a,mu(1))*ones(1,length(Omega)),'--')
loglog(Omega,det_cycle_length(gamma,tau_a,mu(51))*ones(1,length(Omega)),'--')
loglog(Omega,det_cycle_length(gamma,tau_a,mu(end))*ones(1,length(Omega)),'--')
% Generate Fig. 6D:
figure
loglog(mu,T_mean(1,:),'linewidth',15)
hold on
loglog(mu,T_mean(40,:),'linewidth',15)
loglog(mu,T_mean(end,:),'linewidth',15)
loglog(mu,3*tau_a./(Omega(1)*mu))
loglog(mu,3*tau_a./(Omega(40)*mu))
loglog(mu,3*tau_a./(Omega(end)*mu))
loglog(mu,det_cycle_length(gamma,tau_a,mu).*ones(1,length(mu)),'--')
% Generate Fig. 6E:
figure
surf(log10(Omega),log10(mu),log10(T_cv'));
xlim(log10([min(Omega) max(Omega)]))
view(100,20)
% Generate Fig. 6F:
figure
imagesc(log10(mu),log10(Omega),T_cv-1/sqrt(3))
colorbar