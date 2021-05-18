% function void=main_challenge_problem_5b(void)

% Challenge Problem: Stochastic SIR Model -- extension; now estimate little
% r based on the syntheic data
% do this many times - plot curves with deterministic model

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

r_true_value = pars.beta-pars.gamma;


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

% wait to plot the deterministic model after all stochastic runs

%%%% now simulate the stochastic model %%%%%


num_reals = 250; % how many times?

% plot as we go
f1=figure(1); set(gcf,'Position',[240   300   850   400]);
sim_data=struct();
r_estimate=[];
for k =1:num_reals

[this_t_stoch,this_y_stoch] = stochsim_SIR(trange,y0,pars);

% Plot the results
f1=figure(1); subplot(1,2,1);
tmph_stoch_sir(1)=plot(this_t_stoch,this_y_stoch(:,1),'Color',[0.25 0.25 0.25]); hold on;
tmph_stoch_sir(1).Color(4) = 0.1;
tmph_stoch_sir(2)=plot(this_t_stoch,this_y_stoch(:,2),'Color',[0.5 0.5 0.5]); hold on;
tmph_stoch_sir(2).Color(4) = 0.1;
tmph_stoch_sir(3)=plot(this_t_stoch,this_y_stoch(:,3),'Color',[0.75 0.75 0.75]); hold on;
tmph_stoch_sir(3).Color(4) = 0.1;
set(tmph_stoch_sir,'Linewidth',1);

% find time points greater than num days in
num_days=7;
tmpi_num_days=find(this_t_stoch<num_days);

% Fit to incident infections up to 10 days in
[this_p,s]=polyfit(this_t_stoch(tmpi_num_days),log(this_y_stoch(tmpi_num_days,2)),1);


% store the simulated data and estimates
sim_data(k).t_stoch = this_t_stoch;
sim_data(k).S_stoch = this_y_stoch(:,1);
sim_data(k).I_stoch = this_y_stoch(:,2);
sim_data(k).R_stoch = this_y_stoch(:,3);

r_estimate = [r_estimate,this_p(1)];

% p contains the coefficients, s contains the stats
% p
% disp(['this estimate of r = ',num2str(p(1),'%1.4f')]);

end

% Plot results of deterministic SIR model

f1=figure(1); subplot(1,2,1);
% tmph_sir=plot(t,pars.N*y_sir_model); hold on;
tmph_sir(1)=plot(t,pars.N*y_sir_model(:,1),'Color',[0.25 0.25 0.25]); hold on;
tmph_sir(2)=plot(t,pars.N*y_sir_model(:,2),'Color',[0.5 0.5 0.5]); hold on;
tmph_sir(3)=plot(t,pars.N*y_sir_model(:,3),'Color',[0.75 0.75 0.75]); hold on;
set(tmph_sir,'Linewidth',2);
axis([0 60 0 pars.N]);
xlabel('Time (days)');
ylabel('Population fraction');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


% insert legend
% tmpph = [tmph_sir;tmph_si];
tmplh = legend(tmph_sir,'Susceptible','Infectious','Recovered','Location','NorthEast');
legend boxoff;

f1=figure(1); subplot(1,2,2);
h = histogram(r_estimate); hold on;
c = h.BinWidth;
h.BinWidth = 0.01;

f1=figure(1); subplot(1,2,2);
lp(1)=plot(r_true_value*ones(1,10),linspace(0,45,10),'k--','LineWidth',2); hold on;
lp(2)=plot(median(r_estimate)*ones(1,10),linspace(0,45,10),'r','LineWidth',2);
xlabel('estimated r');
ylabel('count');
title([num2str(num_reals),' realizations']);
tmplh = legend(lp,'true value','median estimate','Location','NorthWest');
legend boxoff;

f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';



