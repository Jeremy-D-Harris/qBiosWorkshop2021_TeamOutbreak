function dydt = si_model(t,y,pars)
% SI Model

S=y(1);
I=y(2);

% model equations
dSdt = -pars.beta*S*I+pars.gamma*I; 
dIdt = pars.beta*S*I-pars.gamma*I;

% put in a vector
dydt = [dSdt; dIdt];

