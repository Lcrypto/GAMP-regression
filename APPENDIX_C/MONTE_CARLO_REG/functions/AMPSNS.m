function [beta_means,beta_vars,sigma2,lambdas] = AMPSNS(X,y,maxiter)

% Define prior moments of beta, and alpha
alpha = 10;

% Prior probability for spike and slab prior
lambda = 0.1;
 
% Define prior moments of sigma2 ~ iGamma(a0,b0)
a0 = 0.1;
b0 = 0.1;

% Initialization
[n,p]        =  size(X);
Thresholding =  1.0e-10;
time         =  0;
mu_hat       =  zeros(1,p);       % Initialize vector of parameter estimates at prior mean
tau_hat      =  (1./alpha)*ones(1,p);  % Initialize vector of variances at prior variance
mu_prior     =  ones(1,p);
s            =  zeros(n,1);
pip          =  rand(1,p);


% ====================| ESTIMATION |=================
fprintf('Now you are running AMPSNS')
fprintf('\n')
fprintf('Iteration 0000')
while time < maxiter && norm(mu_hat-mu_prior) > Thresholding
%     clc
    time   =  time + 1;
    if mod(time,500) == 0
        fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%s%4d',8,8,8,8,8,8,8,8,8,8,8,8,8,8,'Iteration ',time)
    end
    mu_prior = mu_hat;
    for jj = 1:1
        % Step 1
        z       = X*mu_hat';
        tau_p   = (X.^2)*tau_hat';
        p_hat   = z - tau_p.*s;
    
        % Set sigma.sq equal to rough estimate of posterior mode
        sigma2  = (2*b0 + (y-z)'*(y-z))/(2*a0 + n);
        
        % Step 2   
        tau_z   = (tau_p.*sigma2)./(tau_p + sigma2);      % Var(z(i) | y,p_hat(k,i),tau_p(k,i))
        z_hat   =  tau_z.*(y./sigma2 + p_hat./tau_p);    % E(z(i) | y,p_hat(k,i),tau_p(k,i))

        s       = (z_hat-p_hat)./tau_p;		
        tau_s   = (1 - tau_z./tau_p) ./ tau_p;
            
        % Step 3
        tau_l   = tau_s'*(X.^2);
        tau_r   = 1./tau_l;
        r_hat   = mu_hat + tau_r.*(s'*X);

        % Step 4
        pip     = 1./( 1 + (1-lambda).*normpdf(zeros(1,p), r_hat, sqrt(tau_r))./(lambda.*normpdf(zeros(1,p), zeros(1,p)-r_hat, sqrt(alpha + tau_l))) );
        nu      = 1./(alpha + tau_l);
        gam     = (r_hat./tau_r).*nu;
        
        mu_hat  = pip.*gam;                               % E(mu(i) | y,r_hat(k,i),tau_r(k,i))
        tau_hat = pip.*(gam.^2 - pip.*gam.^2 + nu);       % Var(mu(i) | y,r_hat(k,i),tau_r(k,i))
    end
    
    lambda = mean(pip);
%     psi    = (1./lambda)*mean( pip.*(gam.^2) + nu);
%     alpha  = 1./psi;
%     alpha(alpha>1e10) = 1e10;
end
fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8)

beta_means   = mu_hat';
beta_vars    = tau_hat';
lambdas      = pip';
