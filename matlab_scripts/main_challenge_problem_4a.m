% function void=main_challenge_problem_1(void)

% Challenge problem: From Epidemic Outbreaks to Endemics

clear all; close all; clc;

% default colors
default_colors = [0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250];

% main data/parameters go here
pars.beta = 0.3; % Transmission rate
pars.gamma = 1/4; % Recovery rate (days^-1)
pars.N = 10000;
pars.I0= 1;
pars.S0= pars.N-pars.I0;
pars.basR0 = pars.beta*pars.S0/pars.gamma/pars.N; % Basic reproduction number

% set up time vector
t_init = 0;
dt = 0.1;
t_end = 365;
pars.t_span = t_init:dt:t_end;


% Run the model
[t,y_sir_model]=ode45(@sir_model,pars.t_span,[pars.S0 pars.I0 0]/pars.N,[],pars);


% Plot the results
f1=figure(1);
tmph_sir=plot(t,y_sir_model,'--'); hold on;
set(tmph_sir,'Linewidth',2);
axis([0 365 0 1]);
xlabel('Time (days)');
ylabel('Population fraction');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

% insert legend
% tmplh = legend('Susceptible','Infectious','Recovered');
tmplh = legend('S (SIR)','Infectious (SIR)','Recovered (SIR)','Location','NorthWest');
legend boxoff;

pause;

% Run the model
[t,y_si_model]=ode45(@si_model,pars.t_span,[pars.S0 pars.I0]/pars.N,[],pars);


% Plot the results
f1=figure(1);
tmph_si(1,1)=plot(t,y_si_model(:,1),'Color',default_colors(1,:)); hold on;
tmph_si(2,1)=plot(t,y_si_model(:,2),'Color',default_colors(2,:));
set(tmph_si,'Linewidth',2);
xlabel('Time (days)');
ylabel('Population fraction');
axis([0 365 0 1]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

% insert legend
tmpph = [tmph_sir;tmph_si];
tmplh = legend('S (SIR)','Infectious (SIR)','Recovered (SIR)','S (SI)','Infectious (SI)','Location','NorthWest');
legend boxoff;


