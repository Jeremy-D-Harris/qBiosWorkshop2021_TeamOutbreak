% function void=main_challenge_problem_2b(void)

% Challenge Problem: Strength-Speed Relationships
% Determine the strength and speed for the following sets of disease parameters

clear all; close all; clc;

% now lets range over the parameters
beta_vary =  [0.5,   1, 0.25, 0.75];
gamma_vary = [0.4, 0.5,  0.5, 0.25];

% main data/parameters go here
% pars.c = 20; % Contacts per unit time (e.g., days) 
% pars.p = 0.025; % Probability of infectious contact
% pars.beta = pars.c*pars.p; % Transmission rate
% pars.gamma = 1/4; % Recovery rate (days^-1)
pars.N = 10000;
pars.I0= 1;
pars.S0= pars.N-pars.I0;
% pars.basR0 = pars.beta*pars.S0/pars.gamma/pars.N; % Basic reproduction number


figure(1); set(gcf,'Position',[240   300   950   400]);

p_array=[];
basR0_array=[];
r_array=[];
for k=1:length(beta_vary)
    
    this_beta = beta_vary(k);
    this_gamma = gamma_vary(k);
    
    pars.beta = this_beta;
    pars.gamma = this_gamma;
    
    this_basR0 = pars.beta*pars.S0/pars.gamma/pars.N;
    basR0_array = [basR0_array; this_basR0];
    
    this_r = pars.beta-pars.gamma;
    r_array = [r_array; this_r];
    
    % set up time vector
    t_init = 0;
    dt = 1; % day intervals
    t_end = 10;
    pars.t_span = t_init:dt:t_end;

    % Run the model over 10 days
    [t,y]=ode45(@sir_model,pars.t_span,[pars.S0 pars.I0 0]/pars.N,[],pars);
    
    % Find the slope
    [this_p,s]=polyfit(t,log(y(:,2)),1);
    
    p_array = [p_array;this_p];
    
    % p contains the coefficients, s contains the stats
    
    figure(1); subplot(1,4,k);
    % Plot the data and overlay the best-fit exponential
    tmph=semilogy(t,y(:,2),'ko','MarkerSize',10); hold on;
    set(tmph,'Linewidth',2); hold on;
    tmph=semilogy(t,exp(this_p(1)*t+this_p(2)),'-','Color',[0.5,0.5,0.5]);
    set(tmph,'linewidth',2);
    
    
    % set up time vector
    t_init = 0;
    dt = 1; % day intervals
    t_end = 30;
    pars.t_span = t_init:dt:t_end;
    
    
    % Use solid points for the future k
    [t,y]=ode45(@sir_model,pars.t_span,[pars.S0 pars.I0 0]/pars.N,[],pars);
    
    % find time points greater than 10 days in
    tmpi=find(t>10);
    
    figure(1); subplot(1,4,k);
    tmph=semilogy(t(tmpi),y(tmpi,2),'ko','MarkerSize',10); hold on;
    set(tmph,'linewidth',2,'markerfacecolor','k');
    xlabel('Time (days)');
    ylabel('Infectious fraction');
    title(['$\mathcal{R}_0 = $',num2str(this_basR0,'%1.1f'),', $r =  $', num2str(this_r,'%1.2f')],'Interpreter','Latex','FontWeight','normal');
    if ismember(k,[2,4])
    axis([0 30 10^-4 1]);
    end
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
end

basR0_actual = basR0_array
% baseR0_estimate = ?

r_actual = r_array
r_estimate = p_array(:,1)