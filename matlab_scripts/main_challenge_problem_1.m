% function void=main_challenge_problem_1(void)

% Challenge Problem: Outbreak Criteria
% Complete the simulation code below for the spread of an infectious
% disease beginning with 1 individual out of 10000 and estimate the value
% of R0. Do you expect the disease to spread or not? Then, change the value
% of p to 0.01 – will the disease spread, why or why not?

clear all; close all; clc;

% main data/parameters go here
pars.c = 20; % Contacts per unit time (e.g., days) 
pars.p = 0.025; % change p to 0.01 % Probability of infectious contact
pars.beta = pars.c*pars.p % Transmission rate
pars.gamma = 1/4; % Recovery rate (days^-1)
pars.N = 10000;
pars.I0= 1;
pars.S0= pars.N-pars.I0;
pars.basR0 = pars.beta*pars.S0/pars.gamma/pars.N % Basic reproduction number

% set up time vector
t_init = 0;
dt = 0.1;
t_end = 100;
pars.t_span = t_init:dt:t_end;


% Run the model
[t,y]=ode45(@sir_model,pars.t_span,[pars.S0 pars.I0 0]/pars.N,[],pars);


% Plot the results
f1=figure(1);
tmph=plot(t,y);
set(tmph,'Linewidth',2);
xlabel('Time (days)');
ylabel('Population fraction');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

% insert legend
tmplh = legend('Susceptible','Infectious','Recovered');
legend boxoff;


