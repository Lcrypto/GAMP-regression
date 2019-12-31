function [beta_means,beta_vars,sigma,lambdas] = SSVS(x,y,nsave,nburn)

T = size(x,1); % time series observations
p = size(x,2); % maximum nubmer of predictors

% ----------------Gibbs related preliminaries
ntot = nsave + nburn;

beta_draws = zeros(nsave,p);
gamma_draws = zeros(nsave,p);
sigma2_draws = zeros(nsave,1);
% ----------------Set priors
% beta ~ N(0,DD), where D = diag(tau_i)
tau_0  = 0.01;
tau_1  = 3;
% gamma_j ~ Bernoulli(1,p_j)
p_j = 0.5*ones(p,1);
% Initialize parameters
gamma = ones(p,1);
h_i = zeros(p,1);

Ftau_0 = (1/tau_0)^2;
Ftau_1 = (1/tau_1)^2;

% Get OLS quanities from the full model (only if degrees of freedom allow this)
if T>p
    beta_OLS = inv(x'*x)'*(x'*y);
    SSE_OLS = (y - x*beta_OLS)'*(y - x*beta_OLS);
    sigma2_OLS = SSE_OLS/(T-(p-1));
end

% Initialize parameters
beta = 0*ones(p,1); %beta_OLS;
sigma2 = 0.1; %sigma2_OLS;


%==========================================================================
%====================| GIBBS ITERATIONS START HERE |=======================
fprintf('Now you are running SSVS')
fprintf('\n')
fprintf('Iteration 0000')
for irep = 1:ntot
    if mod(irep,500) == 0
        fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%s%4d',8,8,8,8,8,8,8,8,8,8,8,8,8,8,'Iteration ',irep)
    end
    
    % 1. Update beta from Normal
    h_i(gamma==1) = Ftau_1;
    h_i(gamma==0) = Ftau_0;
    Delta_beta = inv((x'*x)/sigma2 + diag(h_i));
    miu_beta = Delta_beta*((x'*y)/sigma2);    
    beta = miu_beta + chol(Delta_beta)'*randn(p,1);
    
    % 2. Update restriction indexes of alpha from Bernoulli
    u_i1 = mvnpdf(beta,0,tau_0).*p_j;           
    u_i2 = mvnpdf(beta,0,tau_1).*(1- p_j);
    gst = u_i2./(u_i1 + u_i2);
    gamma = bernoullimrnd(p,gst);       
    
    % 3. Update sigma2 from Inverse Gamma
    c1 = (0.01 + T)/2;
    c2 = 0.5*( 0.01 + (y-x*beta)'*(y-x*beta));
    sigma2 = 1./gamrnd(c1,1./c2);
    
    
    % Save draws
    if irep > nburn
        beta_draws(irep-nburn,:) = beta;
        gamma_draws(irep-nburn,:) = gamma;
        sigma2_draws(irep-nburn,:) = sigma2;
    end
    
end
fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8)
beta_means = squeeze(mean(beta_draws,1))';
beta_vars  = squeeze(var(beta_draws,1))';
lambdas    = squeeze(mean(gamma_draws,1))';
sigma      = squeeze(mean(sigma2_draws,1))';