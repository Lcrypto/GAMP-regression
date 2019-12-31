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
T   = 30;

BETA   = zeros(nMC,T,1,4);
SIGMA  = zeros(nMC,T,4);
ndraws = 1000;          % Number of draws
nburn  = 500;          % Burn in period for MCMC draws

%% Run MC iterations
for iMC = 1:nMC
    disp(['This is iteration ' num2str(iMC)]);
    %% Generate data
    [y,x,beta_sim,sigma_sim] = trend_dgp(T);
    % Estimate
    [T, p] = size(x);
    %x=zscore(x);
    xx=kron(speye(T),ones(1,p));
    xx(xx==1)=x';
    BETA(iMC,:,:,1) = beta_sim; SIGMA(iMC,:,1) = sigma_sim;
    
    % 1) Obtain GAMP
    [beta_SBLTV,~,sigma_SBLTV,al]     = AMPSBL_SV([x, xx],y,1000,1e-10,1e-10);
%     [beta_SBLTV,~,sigma_SBLTV,al]     = AMPSNS_SV([x, xx],y,1000,4,0,0.5,0);
    beta_tvp = zeros(T,p);
    for j=1:p
        beta_tvp(:,j) = beta_SBLTV(j) + beta_SBLTV(p+j:p:end);        
    end
    BETA(iMC,:,:,2) = beta_tvp; SIGMA(iMC,:,2) = sigma_SBLTV;
     
    % 2) Obtain MCMC
    % a) tight prior
    [beta_TVPAR1,sigma_TVPAR1]         = TVP_AR(y,x,ndraws,nburn,0.5);
    BETA(iMC,:,:,3) = squeeze(mean(beta_TVPAR1,3)); SIGMA(iMC,:,3) = squeeze(mean(sigma_TVPAR1,2));
    
    % b) tight prior
    [beta_TVPAR2,sigma_TVPAR2]         = TVP_AR(y,x,ndraws,nburn,10);
    BETA(iMC,:,:,4) = squeeze(mean(beta_TVPAR2,3)); SIGMA(iMC,:,3) = squeeze(mean(sigma_TVPAR2,2));
    
end
save MC1alt_T30.mat
% save(sprintf('%s_%g_%g_%g_%g.mat','MONTE_CARLO',T,p,q,rho),'-mat');
% for i=1:p; subplot(5,4,i); plot([squeeze(mean(BETA(:,:,i,2),1)'), squeeze(mean(BETA(:,:,i,3),1)'), squeeze(mean(BETA(:,:,i,1),1)')],'LineWidth',2);if i==1;legend({'GAMP';'MCMC';'TRUE'});end; end
% figure; plot([mean(BETA(:,:,2),1)', mean(BETA(:,:,3),1)',  beta_sim(:,:)],'LineWidth',2); legend({'GAMP','MCMC','True'});
% figure; plot([mean(SIGMA(:,:,2),1)', mean(SIGMA(:,:,3),1)',mean(SIGMA(:,:,1),1)']); legend({'GAMP','MCMC','True'});
% i=1; plot([squeeze(mean(BETA(:,:,i,2),1)'), squeeze(mean(BETA(:,:,i,3),1)'), squeeze(mean(BETA(:,:,i,4),1)'), squeeze(mean(BETA(:,:,i,1),1)')],'LineWidth',2)

figure;
subplot(3,1,1); plot([squeeze(mean(BETA(1,:,1,2),1)'), squeeze(mean(BETA(1,:,1,1),1)')],'LineWidth',2); grid on; ylim([1 7]);
subplot(3,1,2); plot([squeeze(mean(BETA(1,:,1,3),1)'), squeeze(mean(BETA(1,:,1,1),1)')],'LineWidth',2); grid on; ylim([1 7]);
subplot(3,1,3); plot([squeeze(mean(BETA(1,:,1,4),1)'), squeeze(mean(BETA(1,:,1,1),1)')],'LineWidth',2); grid on; ylim([1 7]);
