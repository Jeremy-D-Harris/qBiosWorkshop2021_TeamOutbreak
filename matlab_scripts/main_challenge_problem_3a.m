% function void=main_challenge_problem_3a(void)

% Challenge Problem: Strength-size in phase plane -- setting up

clear all; close all; clc;

% main data/parameters go here
pars.beta = 0.5; % Transmission rate
pars.gamma = 1/4; % Recovery rate (days^-1)
pars.N = 10000;
pars.I0= 1;
pars.S0= pars.N-pars.I0;


% Modify the initial values
pars.N = 10000;
pars.S0_range=[0.6 0.7 0.8 0.9 0.999];
pars.I0_range=1/pars.N*ones(1,5);
pars.R0_range=1-pars.S0_range-pars.I0_range;
pars.basR0 = pars.beta*pars.S0/pars.gamma/pars.N; % Basic reproduction number


% Run the model
for i=1:length(pars.S0_range)
    
    [t,y]=ode45(@sir_model,[0:0.1:200],[pars.S0_range(i) pars.I0_range(i) pars.R0_range(i)],[],pars);
    tmph=plot(y(:,1),y(:,2),'k-');
    set(tmph,'linewidth',3); hold on;
    tmph=plot(y(end,1),y(end,2),'ko','MarkerSize',10);
    set(tmph,'markerfacecolor',[1 1 1]);
    tmph=plot(y(1,1),y(1,2),'ko','MarkerSize',10);
    set(tmph,'markerfacecolor','k');
end

% Show the excluded regime
patch([1 1 0.75 1],[0 0.25 0.25 0],[0.8 0.8 0.8]);
axis([0 1 0 0.25]);

% % find time points greater than 10 days in
% tmpi=find(t>10);
% tmph=semilogy(t(tmpi),y(tmpi,2),'ko','MarkerSize',10); hold on;
% set(tmph,'linewidth',2,'markerfacecolor','k');
xlabel('Susceptible fraction, \emph{S}','Interpreter','Latex');
ylabel('Infectious fraction, \emph{I}','Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';



