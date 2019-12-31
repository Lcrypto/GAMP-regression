function [beta_means,beta_vars,sigma2,lambdas] = AMPSNS_SV(X,y,maxiter,alpha,upalpha,lambda,uplambda)

% Define prior variance alpha
%alpha = 100;

% Prior probability for spike and slab prior
%lambda = 0.5;

% Components of SMN approximation
pi = [0.0073, .10556, .00002, .04395, .34001, .24566, .2575];
mi = [-10.12999, -3.97281, -8.56686, 2.77786, .61942, 1.79518, -1.08819] - 1.2704;  %% means already adjusted!! %%
sigi = [5.79596, 2.61369, 5.17950, .16735, .64009, .34023, 1.26261];

% Use damping factors
theta_s = 0.9;
theta_mu = 0.9;

% Initialization
[n,p]        =  size(X);
Thresholding =  1.0e-4;
time         =  0;
mu_hat       =  zeros(1,p);       % Initialize vector of parameter estimates at prior mean
tau_hat      =  (1./alpha)*ones(1,p);  % Initialize vector of variances at prior variance
mu_prior     =  ones(1,p);
s            =  zeros(n,1);
pip          =  rand(1,p);
sigma2       =  ones(n,1);

fprintf('Now you are running AMPSNS_SV')
fprintf('\n')
fprintf('Iteration 0000')
% ====================| ESTIMATION |=================
while time < maxiter && norm(mu_hat-mu_prior) > Thresholding
%     clc
    time   =  time + 1;
    if mod(time,100) == 0
        fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%s%4d',8,8,8,8,8,8,8,8,8,8,8,8,8,8,'Iteration ',time)
    end  
    mu_prior = mu_hat;
    for jj = 1:5
        % Step 1
        z       = X*mu_hat';
        tau_p   = (X.^2)*tau_hat';
        p_hat   = z - tau_p.*s;
        
        % Step 2
        tau_z   = (tau_p.*sigma2)./(tau_p + sigma2);     % Var(z(i) | y,p_hat(k,i),tau_p(k,i))
        z_hat   =  tau_z.*(y./sigma2 + p_hat./tau_p);    % E(z(i) | y,p_hat(k,i),tau_p(k,i))

       % s       = (z_hat-p_hat)./tau_p;
        s      =  (1-theta_s)*s + theta_s*( (z_hat-p_hat)./tau_p );
        tau_s   = (1 - tau_z./tau_p) ./ tau_p;
            
        % Step 3
        tau_l   = tau_s'*(X.^2);
        tau_r   = 1./(tau_l + 1e-50);
        r_hat   = mu_hat + tau_r.*(s'*X);

        % Step 4
        l_0     = lnormpdf(zeros(1,p), r_hat, sqrt(tau_r));
        l_1     = lnormpdf(zeros(1,p), zeros(1,p)-r_hat, sqrt(alpha + tau_l));
        pip     = 1./( 1 + ((1-lambda)./lambda).* exp(l_0 - l_1) );
%         pip2     = 1./( 1 + (1-lambda).*(normpdf(zeros(1,p), r_hat, sqrt(tau_r)) + 1e-50)./((lambda.*normpdf(zeros(1,p), zeros(1,p)-r_hat, sqrt(alpha + tau_l))) + 1e-50) );
        nu      = 1./(alpha + tau_l);
        gam     = (r_hat./(tau_r + 1e-50)).*nu;
        
        if sum(isnan(pip))>0; pip(find(isnan(pip))) = 0.5; end
        
        % mu_hat  = pip.*gam;                               % E(mu(i) | y,r_hat(k,i),tau_r(k,i))
        mu_hat = (1-theta_mu)*mu_hat + theta_mu*(pip.*gam);
        tau_hat = pip.*(gam.^2 - pip.*gam.^2 + nu);       % Var(mu(i) | y,r_hat(k,i),tau_r(k,i))
    end
    
    % Set sigma.sq equal to rough estimate of posterior mode
    ly2 = log((y-z).^2 + 1e-6);
    for mix_num = 1:7
        %ls2_mean(:,mix_num) = AMPSBL([ones(n,1), eye(n)],ly2-mi(mix_num),100,1e-10,1e-10);
        ls2_mean(:,mix_num) = ( [mean(ly2-mi(mix_num)); ly2-mi(mix_num)] );
    end
    sigma2 = exp( (sum(repmat(ls2_mean(1,:),n,1).*repmat(pi,n,1),2) + sum(ls2_mean(2:end,:).*repmat(pi,n,1),2))./7 );    
    
    if uplambda == 1   
        lambda = mean(pip);
    end
    if upalpha == 1
        psi    = (1./lambda)*mean( pip.*(gam.^2) + nu);
        alpha  = 1./(psi + 1e-50);
        alpha(alpha>1e10) = 1e10;
        alpha(alpha<10) = 10;
    end
end
fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8)

beta_means   = mu_hat';
beta_vars    = tau_hat';
lambdas      = pip';
