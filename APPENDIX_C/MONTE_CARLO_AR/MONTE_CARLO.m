%% MONTE_CARLO.m
%-------------------------------------------------------------------------------------------
% Written by Dimitris Korobilis
% University of Essex
% This version: February 2017
%-------------------------------------------------------------------------------------------
clear; clc; close all;

addpath('functions')
addpath('MCfunctions')

%% Preliminaries
nMC = 1000;            % Number of Monte Carlo iterations
T   = 500;

BETA   = zeros(nMC,4,3);
ndraws = 1000;          % Number of draws
nburn  = 500;          % Burn in period for MCMC draws

%% Run MC iterations
for iMC = 1:nMC
    disp(['This is iteration ' num2str(iMC)]);
    %% Generate data
    [y,x,beta_sim] = ar_dgp(T);
    % Estimate
    [T, p] = size(x);
    x=zscore(x);
    BETA(iMC,:,1) = beta_sim;
    
    % 1) Obtain GAMP
    [beta_SBL,~,~,~]     = AMPSBL(x,y,1000,1e-10,1e-10);
    BETA(iMC,:,2) = beta_SBL;
     
    % 2) OLS
    BETA(iMC,:,3) = (x'*x)\(x'*y);   
end
save(sprintf('%s_%g.mat','MC_AR',T),'-mat');
