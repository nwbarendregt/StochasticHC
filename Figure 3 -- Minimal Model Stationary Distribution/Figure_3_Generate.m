% Figure_3_Generate.m
% Script used to generate Figure 3 from Barendregt & Thomas, 2021.

clear
% Define model parameters:
alpha = 0.8; beta = 1.3; r = 1;
% Set system-size range and pre-allocate p-value storage:
Omega = logspace(-1,1); p = NaN(1,length(Omega));
% Set sample times and number of desired samples for empirical stationary
% distribution:
dt = 1; t_f = 1e6; t = 0:dt:t_f; N_samples = 1000;
% Define analytical stationary distribution from Eq. (15):
f = @(x,O) (O).^x./(factorial(x)*(exp(O)-1));

for i = 1:length(Omega)
    % Define initial condition (taken once two populations have gone
    % extinct):
    Nint = [ceil(Omega(i));0;0];
    % Simulate realization of minimal model for given system-size:
    [t_sim,N_sim] = GV_gillespie(alpha,beta,r,0,Omega(i),Nint,t_f); N_sim = N_sim(1,:);
    % Pre-allocate state observation vector storage:
    bins = 1:max(N_sim); obsCounts = zeros(1,length(bins));
    % Find observation frequency of each state from realization:
    for j = 1:N_samples
        obsCounts((N_sim(find(t_sim<=t(j+(length(t)-N_samples)),1,'last')))) = obsCounts(N_sim(find(t_sim<=t(j+(length(t)-N_samples)),1,'last')))+1;
    end
    % Calculate expected frequency of each state from analytic
    % distribution:
    expCounts = N_samples*f(bins,Omega(i));
    % Calculate p-value from chi-squared goodness-of-fit test:
    [~,p(i),~] = chi2gof(bins,'ctrs',bins,'frequency',obsCounts,...
        'expected',expCounts);
    % Plot Fig. 1A,B:
    if (i==1) || (i==length(Omega))
        figure
        % Plot histogram of empirical stationary distribution:
        s = stairs(-0.5:(max(N_sim)-0.5),[0 obsCounts/N_samples],'k','linewidth',3);
        hold on
        bot = min(s.YData);
        x = [s.XData(1),repelem(s.XData(2:end),2)];
        y = [repelem(s.YData(1:end-1),2),s.YData(end)];
        fill([x,fliplr(x)],[y,bot*ones(size(y))],[0, 0.4470, 0.7410],'edgecolor',[0, 0.4470, 0.7410])
        % Plot analytical stationary distribution:
        plot([0 bins],[0 f(bins,Omega(i))],'k','linewidth',10)
        xlim([0 max(N_sim(1,:))])
    end
end

% Plot Fig. 1C:
figure
semilogx(Omega,p)
hold on
% Overlay alpha = 0.05 level of significance:
semilogx(Omega,0.05*ones(1,length(Omega)))