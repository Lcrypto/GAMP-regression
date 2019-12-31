function [y_foredraws] = TVP_GK(y,Z,X_fore,constant,nrep,nburn)

% TVP_AR_GK: Time-varying parameters model with mixture innovation
% --------------------------------------------------------------------------
% This code implements the model in Giordani and Kohn (2008).
% ***********************************************************
% The model for inflation (y(t)) is:
%     
%         y(t)  = m(t) +  u(t)
%         m(t) = m(t-1) + K(t)*e(t)
%
%  where u(t)~N(0,h(t)) with stochastic volatility:
%
%         log(h(t)) = log(h(t-1)) + K(t)*eta1(t)
%
%  where eta1(t)~N(0,Q).
% *************************************************************************
%   VERSION: 
%      This version: June 18, 2009
%   COMPATIBILITY: 
%      Written in MATLAB 2008b 64-bit for Windows
%   AUTHOR: 
%      Dimitris Korobilis
%      Department of Economics, University of Strathclyde
%      dikorombilis@yahoo.gr
%
%   Enjoy! (...and extend/modify)
% --------------------------------------------------------------------------

[t,m] = size(Z);
y = y';

% Get OLS quantities
beta_hat = y/Z';
sigma_hat = (y'-Z*beta_hat')'*(y'-Z*beta_hat')./(t-m-1);

%========= PRIORS:
% Use uninformative values
B_OLS = beta_hat'./5;%zeros(m,1);
VB_OLS = sigma_hat*inv(Z'*Z); %#ok<MINV>
h_OLS = .3*ones(m,1);
sigma_OLS = sigma_hat./4;

% Set some hyperparameters here (see page 831, end of section 4.1)
k_Q = 0.01;
k_W = 0.01;

%-------- Now set prior means and variances (_prmean / _prvar)
% B_0 ~ N(B_OLS, 4Var(B_OLS))
B_0_prmean = B_OLS;
B_0_prvar = 4*VB_OLS;
% log(h_0) ~ N(log(h_OLS),I_n)
h_prmean = h_OLS;
h_prvar = 4*eye(m);
% log(sigma_0) ~ N(log(sigma_OLS),I_n)
sigma_prmean = sigma_OLS;
sigma_prvar = 1;

% Note that for IW distribution I keep the _prmean/_prvar notation,
% but these are scale and shape parameters...
% Q ~ IW(k2_Q*size(subsample)*Var(B_OLS),size(subsample))
Q_prmean = ((k_Q).^2)*(1+m)*VB_OLS;
Q_prvar = 1 + m;
% W ~ IW(k2_W*(1+dimension(W))*I_n,(1+dimension(W)))
W_prmean = ((k_W)^2)*2;
W_prvar = 2;

% Parameters of the 7 component mixture approximation to a log(chi^2)
% density:
q_s = [0.00730; 0.10556; 0.00002; 0.04395; 0.34001; 0.24566; 0.25750];
m_s = [-10.12999; -3.97281; -8.56686; 2.77786; 0.61942; 1.79518; -1.08819];
u2_s = [5.79596; 2.61369; 5.17950; 0.16735; 0.64009; 0.34023; 1.26261];

%========= PRIORS ON TRANSISION PROBS (Beta) 
% %1 - Least informative Beta prior
% a_prob = sqrt(t)*0.5;   
% b_prob = sqrt(t)*0.5;

% %2 - Reference Beta prior
% a_prob = 1;   
% b_prob = 1;

%3 - Informative prior (few breaks)
a_prob = 0.1;   
b_prob = 10;

ap_0 = a_prob*ones(2,1);
bp_0 = b_prob*ones(2,1);
% Implied prior "sample size" for state equations
t_0 = (ap_0./(ap_0 + bp_0))*t;

%========= INITIALIZE MATRICES:
% Specify covariance matrices for measurement and state equations
consQ = 0.01;
consW = 0.01;
Wdraw = consW;
Qdraw = consQ*eye(m);
Btdraw = zeros(m,t);
Htdraw = 0.1*ones(m,t-1);
Sigtdraw = 0.1*ones(1,t);
sigtS = ones(t,1);
sigtH = ones(t-1,m);
statedrawS = 5*ones(t,1);
statedrawH = 5*ones(t-1,m);
Zs = ones(t,1);
Zh = kron(ones(t,1),eye(m));
prwH = zeros(numel(q_s),1);
prwS = prwH;

kdraw = 1*ones(t,2);
pdraw = .5*ones(1,2);
kold = kdraw;
kmean = zeros(t,2);
kmax = zeros(t,2);
kvals = ones(2,1);
kvals(1,1) = 0;
kprior = .5*ones(2,1);

% Storage matrices for posteriors and stuff
Bt_draws = zeros(m,t,nrep);
Ht_draws = zeros(m,t-1,nrep);
Sigt_draws = zeros(t,nrep);
Q_draws = zeros(m,m,nrep);
W_draws = zeros(1,nrep);
y_foredraws = zeros(nrep,1);
%----------------------------- END OF PRELIMINARIES ---------------------------

%====================================== START SAMPLING ========================================
%==============================================================================================
tic;
fprintf('Now you are running TVP_GK')
fprintf('\n')
fprintf('Iteration 0000')
for irep = 1:nrep + nburn    % GIBBS iterations starts here
    % Print iterations
    if mod(irep,500) == 0
        fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%s%4d',8,8,8,8,8,8,8,8,8,8,8,8,8,8,'Iteration ',irep)
    end
    
    % I.1: draw K1 index and related probabilities
    %----------------------------------------------------------------------
    ap = ap_0(1,1) + sum(kdraw(:,1));
    bp = bp_0(1,1) + t - sum(kdraw(:,1));

    pdrawa = betarnd(ap,bp);
    pdraw(1,1) = pdrawa;
    kprior(2,1) = pdrawa;
    kprior(1,1) = 1 - kprior(2,1);
    Qchol = chol(Qdraw)';
    
    [kdrawa,~] = gck(y,Z,sqrt(sigtS),kron(ones(t,1),eye(m)),Qchol,kold(:,1),t,zeros(m,1),zeros(m,m),2,kprior,kvals,1,m);
    %kdrawa = 1*ones(t,1);
    kdraw(:,1) = kdrawa;
    kold(:,1) = kdraw(:,1);
    
    % Sample trend inflation (Bt) using Carter and Kohn (1994)
    [Btdraw,~] = carter_kohn_GK(y,Z,sigtS,Qdraw,m,1,t,B_0_prmean,B_0_prvar,kdraw(:,1));
    
    Btemp = Btdraw(:,2:t)' - Btdraw(:,1:t-1)';
    sse_2 = zeros(m,m);
    for i = 1:t-1
        sse_2 = sse_2 + Btemp(i,:)'*Btemp(i,:);
    end
    
    % Qdraw is the variance of the state equation stochastic volatility
    Qinv = inv(sse_2 + Q_prmean);
    Qinvdraw = wish(Qinv,t+Q_prvar);
    Qdraw = inv(Qinvdraw);
   
    % Sample Sigma elements from Normal distribution
    y2 = [];
    for i = 1:t
        ytemps = y(:,i) - Z(i,:)*Btdraw(:,i);
        y2 = [y2  (ytemps.^2)]; %#ok<AGROW>
    end   
    yss = log(y2' + 0.000001);

    %First draw volatilities conditional on sdraw
    vart = zeros(t,1);
    yss1 = zeros(t,1);
    for i = 1:t           
        imix = statedrawS(i,1);
        vart(i,1) = u2_s(imix);
        yss1(i,1) = yss(i,1) - m_s(imix) + 1.2704;
    end
    
    % III.1: draw K2 index and related probabilities
    %----------------------------------------------------------------------
    ap = ap_0(2,1) + sum(kdraw(:,2));
    bp = bp_0(2,1) + t - sum(kdraw(:,2));

    pdrawa = betarnd(ap,bp);
    pdraw(1,2) = pdrawa;
    kprior(2,1) = pdrawa;
    kprior(1,1) = 1 - kprior(2,1);
    Wchol = chol(Wdraw)';
    
    [kdrawa,~] = gck(yss1',2*Zs,vart,ones(t,1),Wchol,kold(:,2),t,0,0,2,kprior,kvals,1,1);
    %kdrawa=1*ones(t,1);
    kdraw(:,2) = kdrawa;
    kold(:,2) = kdraw(:,2);
   
    % Draw the measurement error volatility using CK
    [Sigtdraw,~] = carter_kohn_GK(yss1',2*Zs,vart,Wdraw,1,1,t,sigma_prmean,sigma_prvar,kdraw(:,2));
    
    % Next draw statedraw (chi square approximation mixture component)
    % conditional on Sigtdraw       
    for i = 1:t
        for j = 1:numel(m_s)
            temp1= (1/sqrt(2*pi*u2_s(j)))*exp(-.5*(((yss(i,1) - 2*Sigtdraw(1,i) - m_s(j) + 1.2704)^2)/u2_s(j)));
            prwS(j,1) = q_s(j,1)*temp1;
        end
        prwS = prwS./sum(prwS);
        cprwS = cumsum(prwS);
        trand = rand(1,1);
        if trand < cprwS(1,1); imixS=1;
        elseif trand < cprwS(2,1), imixS=2;
        elseif trand < cprwS(3,1), imixS=3;
        elseif trand < cprwS(4,1), imixS=4;
        elseif trand < cprwS(5,1), imixS=5;
        elseif trand < cprwS(6,1), imixS=6;
        else imixS=7; 
        end
        statedrawS(i,1)=imixS;
    end
        
    sigtS = exp(Sigtdraw');

    Sigttemp = Sigtdraw(:,2:t)' - Sigtdraw(:,1:t-1)';

    sse_2 = 0;
    for i = 1:t-1
        sse_2 = sse_2 + Sigttemp(i,:)'*Sigttemp(i,:);
    end
    
    % Draw the variance of the meaasurement error stochastic volatility
    Winv = inv(sse_2 + W_prmean);
    Winvdraw = wish(Winv,t+W_prvar);
    Wdraw = inv(Winvdraw);
    
    %----------------------------SAVE AFTER-BURN-IN DRAWS AND IMPULSE RESPONSES -----------------
    if irep > nburn
        Bt_draws(:,:,irep-nburn) = Btdraw;
        Ht_draws(:,:,irep-nburn) = Htdraw;
        Sigt_draws(:,irep-nburn) = Sigtdraw';
        Q_draws(:,:,irep-nburn) = Qdraw;
        W_draws(:,irep-nburn) = Wdraw;
        
        %==========Forecasting
        betas_fore = [Btdraw(:,t)' ; zeros(1,m)];
        sigma2_fore = [sigtS(t); zeros(1,1)];                   
        % Update the parameters
        betas_fore(2,:) = betas_fore(1,:) + randn(1,m)*chol(Qdraw);
        % We do not want to forecast with explosive
        % parameter draws. Here I check only the AR(1)
        % coefficient and if it is explosive I set it back
        % to the previous regime (instead of drawing again
        % until stationary which can be inefficient computationally)
        if sum(abs(betas_fore(2,1*constant + 1:end))) >= 1.000                        
            betas_fore(2,1*constant + 1) = betas_fore(1,1*constant + 1);
        end
        log_sigma_tplusone = log(sigma2_fore(1,:)) + sqrt(Wdraw)*randn;
        sigma2_fore(2,:) = exp(log_sigma_tplusone);    
        
        y_hat = X_fore*betas_fore(2,:)' + sqrt(sigma2_fore(2))*randn;
        y_foredraws(irep-nburn,:) = y_hat;  
    end  % End saving after burn-in draws
end %END main Gibbs loop (for irep = 1:nrep+nburn)
fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8)
toc;
