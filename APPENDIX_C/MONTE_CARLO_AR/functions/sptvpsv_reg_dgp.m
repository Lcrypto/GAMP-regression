function [y,x,s,theta_t,beta_t,sigma_t] = sptvpsv_reg_dgp(T,p,param)

% **************************************************************************************************************************
% SPTVP_REG_DGP  Code that generates data from a univariate sparse, time-varying parameter (tvp) regression with
% constant variance. 
% The model is of the following form
%
%                   y[t]  =  x[t] beta[t] + e[t]
%                beta[t]  =  s[t] theta[t]
%               theta[t]  =  mu + 0.99I x (theta[t-1] - mu) + u[t]
%
% where e[t] ~ N(0, sigma), u[t] ~ N(0,(1/T)^2 x I), mu is the unconditional mean of the AR process for theta[t], s[t] is a 
% (px1) vector of 0/1 values (sparsity indicator), and we also assume the initial condition theta[0] to be fixed and
% known when generating this process.
% **************************************************************************************************************************
% INPUTS:
%  T        Time series observations
%  p        Number of predictors   
%  param    Structure array that specified all other tuning parameters
%      param.rho     Correlation coefficient for predictors x
%      param.q       Number of predictors that have non-zero coefficients (i.e. they are important/"significant" for y
%      param.s       Generate the (Txp) variable s, which indexes which coefficients are zero (or not) in what periods
%      param.sigma0  Initial condition for the process for log(sigma[t])
%      param.theta0  Initial condition for the process for theta[t]
%      param.mu      Unconditional mean for the process for theta[t]
%
% OUTPUTS:
%  y        Generated time series following the sparse tvp regression
%  x        Generated right-hand side predictor variables
%  theta_t  Generated unrestricted coefficients theta[t]
%  beta_t   Generated restricted (sparse) coefficients beta[t]
%  sigma_t    Regression variance
%
% **************************************************************************************************************************
% Written by Dimitris Korobilis on 22/04/2017
% University of Essex
% **************************************************************************************************************************

%% Check for INPUT arguments
if nargin == 0
    T   = 200;          % Time series observations
    p   = 8;            % Number of predictors
    rho = 0;            % Correlation between predictors is rho^|i-j| 
    q   = round(.5*p);  % Percentage of non-zero initial parameters
    
    s   = zeros(T,p);
    s(:,[1 3 5 7]) = 1;
%     s(1:round(T/2),[1 3 5 7 9 10]) = 1;
%     s(round(T/2)+1:end,[1 2 4 5 7 10]) = 1;    
%     s(1:round(T/2),randperm(p,q)) = 1;
%     s(round(T/2)+1:end,randperm(p,q)) = 1;
    
    V0 = log(.2);  % Regression variance
    theta0 = 2.*rand(1,p)+0.5;%[1.5, 0, -1, 0, 1.2,  0, 3.2, 0];  % Initial regression coefficients
    mu = s(1,:).*theta0;%[1.5, 0, -1, 0, 1.2,  0, 3.2, 0];       % Long-run (unconditional) mean of regression coefficients
else
    rho = param.rho;
    s = param.s;
    V0 = log(param.sigma0); theta0 = param.theta0; mu = param.mu;
end

%% Generate sparse time-varying parameter model
% 1. First generate (possibly) correlated predictors
corr_x = zeros(p,p);          % Correlation matrix for predictors
for i = 1:p                   % do the lazy version, using for loops
    for j = 1:p           
        corr_x(i,j) = rho^(abs(i-j));
    end
end
x = randn(T,p);%*chol(corr_x);  % Generate RHS predictors

% 2. Generate unrestricted regression coefficients theta_t
theta_t = zeros(T+100,p);
for t = 1:T+100
    if t == 1
        theta_t(t,:) = mu +  0.99*(theta0 - mu) + (1/(T^(3/4)))*randn(1,p);
    else
        theta_t(t,:) = mu +  0.99*(theta_t(t-1,:) - mu) + (1/(T^(3/4)))*randn(1,p);
    end
end
theta_t = theta_t(101:end,:);

% 3. Define sparse (restricted) regression coefficients beta_t
beta_t = s.*theta_t;


% 4. Generate unrestricted regression coefficients sigma_t
V_t = zeros(T+100,1);
for t = 1:T+100
    if t == 1
        V_t(t,:) = V0 +  0.98*(2*V0 - V0) + (1/(T^(2/4)))*randn;
    else
        V_t(t,:) = V0 +  0.98*(V_t(t-1,:) - V0) + (1/(T^(2/4)))*randn;
    end
end
sigma_t = exp(V_t(101:end,:));

% 5. Generate dependent variable y
y = sum(x.*beta_t,2) + sqrt(sigma_t).*randn(T,1);

