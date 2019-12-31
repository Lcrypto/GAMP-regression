function [beta_means,beta_vars,sigma2,lambdas] = AMPSBL_SV(X,y,maxiter,a,b)

% Define prior moments of alpha ~ Gamma(a,b)
%a = 1e-1;%1-0.2*size(X,1)/size(X,2);
%b = 1e-1;

% Components of SMN approximation
pi = [0.0073, .10556, .00002, .04395, .34001, .24566, .2575];
mi = [-10.12999, -3.97281, -8.56686, 2.77786, .61942, 1.79518, -1.08819] - 1.2704;  %% means already adjusted!! %%
sigi = [5.79596, 2.61369, 5.17950, .16735, .64009, .34023, 1.26261];

% Initialization
[n,p]        =  size(X);
Thresholding =  1.0e-10;
time         =  0;
mu_hat       =  zeros(p,1);  % in high dimensional case. replaced as At*y since A_t*A is highly computing cost
tau_hat      =  ones(p,1);
mu_prior     =  ones(p,1);
s            =  zeros(n,1); 
alpha        =  ones(p,1);
sigma2       =  zeros(n,1);

fprintf('Now you are running AMPSBL_SV')
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
        ly2 = log((y-z).^2 + 1e-10);

        for mix_num = 1:7
            ls2_mean(:,mix_num) = ( [mean(ly2-mi(mix_num)); ly2-mi(mix_num)] );
        end        
        sigma2 = exp(sum(ls2_mean(1,:).*repmat(pi,n,1),2)) + exp(sum(ls2_mean(2:end,:).*repmat(pi,n,1),2))./7;
        
        % Step 2
        tau_z  =  (tau_p.*sigma2)./(tau_p + sigma2);       % Var(z(i) | y,p_hat(k,i),tau_p(k,i))
        z_hat  =  tau_z.*(y./sigma2 + p_hat./tau_p);      % E(z(i) | y,p_hat(k,i),tau_p(k,i))
        
        s      =  (z_hat-p_hat)./tau_p;
        tau_s  =  (1 - tau_z./tau_p)./tau_p;
    
        % Step 3 (variable update part)
        tau_l  =  (X.^2)'*tau_s;
        tau_r  =  1./(tau_l + 1e-50);
        r_hat  =  mu_hat + tau_r.*(X'*s);
    
        % Step 4
        mu_hat = (r_hat.*tau_l)./(alpha + tau_l);
        tau_hat = 1./(alpha + tau_l);
    end
    
    % En (parameters update)
    alpha = (a + .5)./(0.5*(mu_hat.^2) + b); %prior variance for the coeff beta
    alpha(alpha>1e15) = 1e15;    
end

fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8)

beta_means   = mu_hat;
beta_vars    = tau_hat;
lambdas      = alpha;