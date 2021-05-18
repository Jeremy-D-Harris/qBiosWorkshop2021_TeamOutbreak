function dydt = sir_model(t,y,pars)
% SIR Model

S=y(1);
I=y(2);

% model equations
dSdt = -pars.beta*S*I;
dIdt = pars.beta*S*I-pars.gamma*I; 
dRdt = pars.gamma*I;

% put in a vector
dydt = [dSdt; dIdt; dRdt];