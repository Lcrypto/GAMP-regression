function [y,x,s,theta_t,beta_t,sigma] = sptvp_reg_dgp(T,p,param)

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
%      param.sigma   Variance of y, i.e. variance of the regression model
%      param.theta0  Initial condition for the process for theta[t]
%      param.mu      Unconditional mean for the process for theta[t]
%
% OUTPUTS:
%  y        Generated time series following the sparse tvp regression
%  x        Generated right-hand side predictor variables
%  theta_t  Generated unrestricted coefficients theta[t]
%  beta_t   Generated restricted (sparse) coefficients beta[t]
%  sigma    Regression variance
%
% **************************************************************************************************************************
% Written by Dimitris Korobilis on 22/04/2017
% University of Essex
% **************************************************************************************************************************

%% Check for INPUT arguments
if nargin == 0
    T   = 200;            % Time series observations
    p   = 10;            % Number of predictors
    rho = 0;           % Correlation between predictors is rho^|i-j| 
    q   = round(.5*p);  % Percentage of non-zero initial parameters
    
    s   = zeros(T,p);
    s(1:round(T/2),randperm(p,q)) = 1;
    s(round(T/2)+1:end,randperm(p,q)) = 1;
    
    sigma = 1;  % Regression variance
    theta0 = 4*rand(1,p);  % Initial regression coefficients
    mu = 0.5*theta0;       % Long-run (unconditional) mean of regression coefficients
else
    rho = param.rho; q = param.q;
    s = param.s;
    sigma = param.sigma; theta0 = param.theta0; mu = param.mu;
end

%% Generate sparse time-varying parameter model
% 1. First generate (possibly) correlated predictors
corr_x = zeros(p,p);          % Correlation matrix for predictors
for i = 1:p                   % do the lazy version, using for loops
    for j = 1:p           
        corr_x(i,j) = rho^(abs(i-j));
    end
end
x = randn(T,p)*chol(corr_x);  % Generate RHS predictors

% 2. Generate unrestricted regression coefficients theta_t
theta_t = zeros(T+50,p);
for t = 1:T+50
    if t == 1
        theta_t(t,:) = mu +  0.99*(theta0 - mu) + (1/(T^(3/4)))*randn(1,p);
    else
        theta_t(t,:) = mu +  0.99*(theta_t(t-1,:) - mu) + (1/(T^(3/4)))*randn(1,p);
    end
end
theta_t = theta_t(51:end,:);

% 3. Define sparse (restricted) regression coefficients beta_t
beta_t = s.*theta_t;

% 4. Generate dependent variable y
y = sum(x.*beta_t,2) + randn(T,1)*sqrt(sigma);

