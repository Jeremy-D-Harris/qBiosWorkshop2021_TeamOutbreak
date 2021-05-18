% function void=main_challenge_problem_2a(void)

% Challenge Problem: Strength-Speed Relationships -- setting up
% Determine the strength and speed for the following sets of disease parameters

% beta =  0.5,   1, 0.25, 0.75
% gamma = 0.4, 0.5,  0.5, 0.25


clear all; close all; clc;

% main data/parameters go here
pars.c = 20; % Contacts per unit time (e.g., days) 
pars.p = 0.025; % Probability of infectious contact
pars.beta = pars.c*pars.p; % Transmission rate
pars.gamma = 1/4; % Recovery rate (days^-1)
pars.N = 10000;
pars.I0= 1;
pars.S0= pars.N-pars.I0;
pars.basR0 = pars.beta*pars.S0/pars.gamma/pars.N; % Basic reproduction number

% set up time vector
t_init = 0;
dt = 1; % day intervals
t_end = 10;
pars.t_span = t_init:dt:t_end;


% Run the model over 10 days 

[t,y]=ode45(@sir_model,pars.t_span,[pars.S0 pars.I0 0]/pars.N,[],pars);

% Find the slope
[p,s]=polyfit(t,log(y(:,2)),1);

% p contains the coefficients, s contains the stats
p

% Plot the data and overlay the best-fit exponential 
tmph=semilogy(t,y(:,2),'ko','MarkerSize',10); hold on;
set(tmph,'Linewidth',2); hold on;
tmph=semilogy(t,exp(p(1)*t+p(2)),'-','Color',[0.5,0.5,0.5]); 
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
tmph=semilogy(t(tmpi),y(tmpi,2),'ko','MarkerSize',10); hold on;
set(tmph,'linewidth',2,'markerfacecolor','k');
xlabel('Time (days)');
ylabel('Population fraction');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';



