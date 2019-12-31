%% MONTE_CARLO.m Monte Carlo exercise to replicate the results in Korobilis and Ribeiro (2017)
%--------------------------------------------------------------------------------------------
%  Generate data from regression model and compare the following estimators
%  1) AMP - LASSO
%  2) AMP - Spike and Slab
%  3) Gibbs -  LASSO
%  4) Gibbs - Spike and Slab
%-------------------------------------------------------------------------------------------
% Written by Dimitris Korobilis
% University of Essex
% This version: February 2017
%-------------------------------------------------------------------------------------------

clear; close all; clc;

addpath('functions')

%% Preliminaries
nMC = 500;            % Number of Monte Carlo iterations
T   = 50;             % Number of observations to generate
p   = 100;            % Number of predictors to generate
q   = round(0.01*p);  % Number of important predictors
rho = 0.3;            % Correlation among predictors
SNR = 1;              % Signal-to-Noise Ratio
%% Start MC loop
BETA   = zeros(nMC,p,6);
Etimes = zeros(nMC,4);
for iMC = 1:nMC
    disp(['Now running Monte Carlo Iteration ' num2str(iMC)]);
    %% -----------------------------Generate artificial data--------------------------------
    corr_x = zeros(p,p);  % Correlation matrix for predictors
    for i = 1:p           % do the lazy version, using for loops
        for j = 1:p
            corr_x(i,j) = rho^(abs(i-j));
        end
    end
    x = randn(T,p)*chol(corr_x);              % Generate RHS predictors
    beta_sim = [4*(2*rand(q,1)-1); zeros(p-q,1)]; % Regression coefficients
    sigma_sim = 1;%(beta_sim'*corr_x*beta_sim)./SNR;        % Regression variance
    %now you are ready to simulate y
    y = x*beta_sim + sqrt(sigma_sim)*randn(T,1);
    %----------------------------------------------------------------------------------------
    
    
    %% Begin estimation using the various methods    
    % First use uncorrelated predictors, if we have less predictors than observations
    if p<T
        r = cholcov(cov(x)); %r = sqrtm(cov(x));
        xu=x/r;
    else        
        xu = zscore(x);   r = eye(p);
    end

    % AMP Sparse Bayesian Learning (SBL)
    tic;
    [beta_SBL,~,~,alpha] = AMPSBL(xu,y,1000);
    t_sbl = toc;
    beta_SBL = ((r'*r)\r')*beta_SBL;
    
    % AMP Spike and Slab (SNS)
    tic;
    [beta_SNS,~,~,pip] = AMPSNS(xu,y,1000);
    t_sns = toc;
    beta_SNS = ((r'*r)\r')*beta_SNS;

    % Gibbs sampler LASSO
    tic;
    [beta_LASSO,~,~,lambda2] = LASSO(xu,y,2000,1000);
    t_lasso = toc;
    beta_LASSO = ((r'*r)\r')*beta_LASSO;
    
    % Gibbs sampler Spike and Slab
    tic;    
    [beta_SSVS,~,~,lambda2] = SSVS(xu,y,2000,1000);
    t_ssvs = toc;
    beta_SSVS = ((r'*r)\r')*beta_SSVS;
    
    for i = 1:p
        beta_OLS(i,:) = (xu(:,i)'*xu(:,i))\(xu(:,i)'*y);
    end
    beta_OLS = ((r'*r)\r')*beta_OLS;
        
    beta = [beta_sim, beta_OLS, beta_SBL, beta_SNS, beta_LASSO, beta_SSVS];
    BETA(iMC,:,:) = beta;
    Etimes(iMC,:) = [t_sbl,t_sns,t_lasso,t_ssvs];
end

save(sprintf('%s_%g_%g_%g.mat','MONTE_CARLO',T,p,rho),'-mat');
