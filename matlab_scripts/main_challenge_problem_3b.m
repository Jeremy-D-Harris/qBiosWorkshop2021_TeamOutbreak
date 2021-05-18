% function void=main_challenge_problem_3b(void)

% Challenge Problem: Strength-size in phase plane
% First, fix the value of the recovery rate to γ = 0.25, then modulate the
% transmission rate from β = 0.25 to β = 2.5. Assume that initially, 80% of
% the population is susceptible and 20% is recovered. In this case, does
% the disease always spread? Why or why not? In addition, find the outbreak
% size and show how it relates to Reff – the effective reproduction number
% given the initial susceptible population.

clear all; close all; clc;

% main data/parameters go here
pars.beta = 0.5;
pars.N=10000;
pars.beta_range=[0.25 0.3 0.35 0.5 1 1.5 2.5];  % Transmission rate
pars.gamma = 1/4; % Recovery rate (days^-1)
pars.S0=0.8;
pars.I0=1/pars.N;
pars.R0=1-pars.S0-pars.I0;

figure(1); set(gcf,'Position',[240   300   950   400]);

% Run the model - varying over beta range
for i=1:length(pars.beta_range)
    
    pars.beta=pars.beta_range(i);
    pars.effR0(i) = pars.beta/pars.gamma*pars.S0; 
    pars.r(i)=pars.gamma*(pars.effR0(i)-1); 
    [t,y]=ode45(@sir_model,[0:0.1:200],[pars.S0 pars.I0 pars.R0],[],pars); 
    pars.size(i)=y(end,3)-pars.R0;
    
    figure(1); subplot(1,2,1);
    tmph=plot(y(:,1),y(:,2),'k-');
    set(tmph,'linewidth',2,'color',[0.75 0.75 0.75] - [0.1*i 0.1*i 0.1*i]); hold on
    tmph=plot(y(end,1),y(end,2),'ko','MarkerSize',10);
    set(tmph,'markerfacecolor',[1 1 1]);
    tmph=plot(y(1,1),y(1,2),'ko','MarkerSize',10); hold on;
    set(tmph,'markerfacecolor','k');
    
end


% Show the excluded regime
patch([1 1 0.5 1],[0 0.5 0.5 0],[0.8 0.8 0.8]);
axis([0 1 0 0.5]);
xlabel('Susceptible fraction, \emph{S}','Interpreter','Latex');
ylabel('Infectious fraction, \emph{I}','Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


% pause;

% main data goes here 
pars.N=10000; 
pars.beta_range=[0.01:0.01:2.5]; 
pars.gamma = 1/4;
pars.S0=0.8; 
pars.I0=1/pars.N; 
pars.R0=1-pars.S0-pars.I0;

% Run the model
for i=1:length(pars.beta_range)
pars.beta=pars.beta_range(i);
pars.effR0(i) = pars.beta/pars.gamma*pars.S0; 
pars.r(i)=pars.gamma*(pars.effR0(i)-1); 
[t,y]=ode45(@sir_model,[0:0.1:200],[pars.S0 pars.I0 pars.R0],[],pars); 
pars.size(i)=y(end,3)-pars.R0;
end

% show effective reproduction numbers
figure(1); subplot(1,2,2);
tmph=plot(pars.beta_range/pars.gamma,pars.size,'k-'); 
set(tmph,'linewidth',2); hold on;
tmph=plot([1/pars.S0 1/pars.S0],[0 1],'k--'); 
set(tmph,'linewidth',2);
set(gca,'xtick',[0 1/pars.S0 2 4 6 8 10]);
xlabel('Transmission/recover, $\beta/\gamma$','Interpreter','Latex');
ylabel('Outbreak Size, $R_{\infty} - R_{0}$','Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

