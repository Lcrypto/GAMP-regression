function [beta_means,beta_vars,sigma2,lambdas] = AMPSBL(X,y,maxiter)

% Define prior moments of alpha ~ Gamma(a,b)
a = 1e-10;%1-0.2*size(X,1)/size(X,2);
b = 1e-10;%1;

% Define prior moments of sigma2 ~ iGamma(a0,b0)
a0 = 0.1;
b0 = 0.1;

% Initialization
[n,p]        =  size(X);
Thresholding =  1.0e-10;
time         =  0;
mu_hat       =  zeros(p,1);  % in high dimensional case. replaced as At*y since A_t*A is highly computing cost
tau_hat      =  ones(p,1);
mu_prior     =  ones(p,1);
s            =  zeros(n,1); 
alpha        =  ones(p,1);
fprintf('Now you are running AMPSBL')
fprintf('\n')
fprintf('Iteration 0000')
% iterative process
while time < maxiter && norm(mu_hat-mu_prior) > Thresholding
%     clc
    time   =  time + 1;
    if mod(time,500) == 0
        fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%s%4d',8,8,8,8,8,8,8,8,8,8,8,8,8,8,'Iteration ',time)
    end   
    mu_prior = mu_hat;
    for jj=1:1
        % Step 1 (factor update part)
        z      =  X*mu_hat;
        tau_p  =  (X.^2)*tau_hat;
        p_hat  =  z - tau_p.*s;
    
        % Set sigma.sq equal to rough estimate of posterior mode   
        sigma2 = (2*b0 + (y-z)'*(y-z))/(2*a0 + length(y));
        
        % Step 2
        tau_z  =  (tau_p*sigma2)./(tau_p + sigma2);       % Var(z(i) | y,p_hat(k,i),tau_p(k,i))
        z_hat  =  tau_z.*(y./sigma2 + p_hat./tau_p);      % E(z(i) | y,p_hat(k,i),tau_p(k,i))
        
        s      =  (z_hat-p_hat)./tau_p;
        tau_s  =  (1 - tau_z./tau_p)./tau_p;
    
        % Step 3 (variable update part)
        tau_l  =  (X.^2)'*tau_s;
        tau_r  =  1./tau_l;
        r_hat  =  mu_hat + tau_r.*(X'*s);
    
        % Step 4
        mu_hat = (r_hat.*tau_l)./(alpha + tau_l);
        tau_hat = 1./(alpha + tau_l);
    end
    
    % En (parameters update)
    alpha = (a -.5)./(0.5*(mu_hat.^2) + b); %prior variance for the coeff beta
    alpha(alpha>1e10) = 1e10;    
end
fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8)

beta_means   = mu_hat;
beta_vars    = tau_hat;
lambdas      = alpha;