function [t,y] = stochsim_SIR(trange,y0,pars)

% Simulates an SIR model via the Gillespie algorithm % from t0 to tf in
% trange given initial conditions in y0 = [S0 I0 R0] and parameters in
% pars. Returns time and values

% be careful with beta and gamma
beta_I = pars.beta;
gamma_I = pars.gamma;
N = pars.N;

% Conditions
t0=trange(1);
tf=trange(2);
t(1)=t0;
y(1,:)=y0;
tcur=t0;
ycur=y0;
ind=1;

% Model
while (tcur<tf)
    
    
    % y = [S I R];
    S_curr = ycur(1);
    I_curr = ycur(2);
    R_curr = ycur(3);
    
    % Check to see if there is an infection
    if (I_curr==0)
        ind=ind+1;
        t(ind)=tf;
        y(ind,:)=ycur;
%         disp('no more infections \n');
        break;
    end
    
    % Rates
    infrate = beta_I*S_curr*I_curr/N; % infection rate
    recrate = gamma_I*I_curr; % recovery rate
    totrate = infrate + recrate; % total rate
    
    % get next time
    dt = -1/totrate*log(rand);
    tcur=tcur+dt;
    
    % Event type - which event happens
    if (rand<(infrate/totrate)) % infection
        
        % update variables
        S_curr = S_curr-1; % one less susceptible
        I_curr = I_curr+1; % one more infected individual
        
        ycur = [S_curr,I_curr,R_curr];
        
    else % recovery happens
        
        % update variables
        I_curr = I_curr-1; % one less infected individual
        R_curr = R_curr+1; % one more recovered person
        
        ycur = [S_curr,I_curr,R_curr];
        
    end
    
    ind=ind+1; % indexing
    t(ind)=tcur; % time series
    y(ind,:)=ycur; % states over time
    
end

