% function void=main_challenge_problem_5a(void)

% Challenge Problem: Stochastic SIR Model -- set up

clear all; close all; clc;

% default colors
default_colors = [0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250];

%%%%% first simulate the deterministic model %%%%%

% main data/parameters go here
pars.c = 20; % Contacts per unit time (e.g., days) 
pars.p = 0.025; % Probability of infectious contact
pars.beta = pars.c*pars.p; % = 0.5 - Transmission rate
pars.gamma = 1/4; % Recovery rate (days^-1)
pars.N = 1000;
pars.I0= 10; % start with 10 infections
pars.S0= pars.N-pars.I0;
pars.basR0 = pars.beta*pars.S0/pars.gamma/pars.N; % Basic reproduction number

% set up time vector
t_init = 0;
dt = 0.1;
t_end = 60;
pars.t_span = t_init:dt:t_end;
trange = [t_init t_end];

% initial condition
y0 = [pars.S0 pars.I0 0]; % numbers

% First run the deterministic SIR model
[t,y_sir_model]=ode45(@sir_model,pars.t_span,y0/pars.N,[],pars);


% Plot results of deterministic SIR model
f1=figure(1);
tmph_sir=plot(t,pars.N*y_sir_model); hold on;
set(tmph_sir,'Linewidth',2);
axis([0 60 0 pars.N]);
xlabel('Time (days)');
ylabel('Population fraction');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


%% now simulate the stochastic model %%%%%
[t_stoch,y_stoch] = stochsim_SIR(trange,y0,pars);

% Plot the results
f1=figure(1);
tmph_stoch_sir(1)=plot(t_stoch,y_stoch(:,1),'--','Color',default_colors(1,:)); hold on;
tmph_stoch_sir(2)=plot(t_stoch,y_stoch(:,2),'--','Color',default_colors(2,:)); hold on;
tmph_stoch_sir(3)=plot(t_stoch,y_stoch(:,3),'--','Color',default_colors(3,:)); hold on;
set(tmph_stoch_sir,'Linewidth',2);

% find time points greater than 10 days in
num_days=7;
tmpi_num_days=find(t_stoch<num_days);

% Fit to incident infections up to 10 days in
[p,s]=polyfit(t_stoch(tmpi_num_days),log(y_stoch(tmpi_num_days,2)),1);

% p contains the coefficients, s contains the stats
% p
disp(['this estimate of r = ',num2str(p(1),'%1.4f')]);

% Plot the data and overlay the best-fit exponential 
f2 = figure(2);
tmph=semilogy(t_stoch(tmpi_num_days),y_stoch(tmpi_num_days,2),'ko','MarkerSize',10); hold on;
set(tmph,'Linewidth',1); hold on;
tmph=semilogy(t_stoch,exp(p(1)*t_stoch+p(2)),'-','Color',[0.5,0.5,0.5]); 
set(tmph,'linewidth',2);


% find time points greater than 10 days in
tmpi=find(t_stoch>num_days);
tmph=semilogy(t_stoch(tmpi),y_stoch(tmpi,2),'ko','MarkerSize',10); hold on;
set(tmph,'linewidth',2,'markerfacecolor',[0 0 0]);
xlabel('Time (days)');
ylabel('Infected population');
axis([0 21 10 pars.N]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';





