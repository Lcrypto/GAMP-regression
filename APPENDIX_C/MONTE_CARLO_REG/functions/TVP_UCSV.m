function [Bt_draws,Sigt_draws] = TVP_UCSV(y,nrep,nburn)
% TVP_AR_SW: Unobserved components model with stochastic volatility
% --------------------------------------------------------------------------
% This code implements the model in Stock and Watson (2007).
% ***********************************************************
% The model for inflation (y(t)) is:
%     
%         y(t)  = m(t) +  u(t)
%         m(t) = m(t-1) + e(t)
%
%  where u(t)~N(0,h(t)), e(t)~N(0,w(t)) with stochastic volatilities:
%
%         log(h(t)) = log(h(t-1)) + eta1(t)
%         log(w(t)) = log(w(t-1)) + eta2(t)
%
%  where eta1(t)~N(0,Q) and eta2~N(0,W).
% *************************************************************************
%   VERSION: 
%      This version: May 20, 2010
%   COMPATIBILITY: 
%      Written in MATLAB 2008b 64-bit for Windows
%   AUTHOR: 
%      Dimitris Korobilis
%      Department of Economics, University of Strathclyde
%      dikorombilis@yahoo.gr
% --------------------------------------------------------------------------

m=1; p=1;
t = size(y,1);
Z = ones(t,1);
y=y';
%========= PRIORS:
%this code was initially set up to use OLS on a training sample to calibrate prior
%but let us just use relatively noninformative values
B_OLS = zeros(m,1);
VB_OLS = eye(m);
h_OLS = ones(m,1);
sigma_OLS = 1;

% Set some hyperparameters here (see page 831, end of section 4.1)
k_Q = 0.1;
k_W = 0.1;

%-------- Now set prior means and variances (_prmean / _prvar)
% B_0 ~ N(B_OLS, 4Var(B_OLS))
B_0_prmean = B_OLS;
B_0_prvar = 4*VB_OLS;
% log(h_0) ~ N(log(h_OLS),I_n)
h_prmean = h_OLS;
h_prvar = 4*eye(m);
% log(sigma_0) ~ N(log(sigma_OLS),I_n)
sigma_prmean = sigma_OLS;
sigma_prvar = 4;

% Note that for IW distribution I keep the _prmean/_prvar notation,
% but these are scale and shape parameters...
% Q ~ IW(k2_Q*size(subsample)*Var(B_OLS),size(subsample))
Q_prmean = ((k_Q).^2)*(1+m)*VB_OLS;
Q_prvar = 1 + m;
% W ~ IW(k2_W*(1+dimension(W))*I_n,(1+dimension(W)))
W_prmean = ((k_W)^2)*2*eye(p);
W_prvar = 2;

% Parameters of the 7 component mixture approximation to a log(chi^2)
% density:
q_s = [0.00730; 0.10556; 0.00002; 0.04395; 0.34001; 0.24566; 0.25750];
m_s = [-10.12999; -3.97281; -8.56686; 2.77786; 0.61942; 1.79518; -1.08819];
u2_s = [5.79596; 2.61369; 5.17950; 0.16735; 0.64009; 0.34023; 1.26261];

%========= INITIALIZE MATRICES:
% Specify covariance matrices for measurement and state equations
consQ = 0.000001;
consW = 0.0001;
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

% Storage matrices for posteriors and stuff
Bt_draws = zeros(t,nrep);
Sigt_draws = zeros(t,nrep);
%----------------------------- END OF PRELIMINARIES ---------------------------

%====================================== START SAMPLING ========================================
%==============================================================================================
tic; % This is just a timer
fprintf('Now you are running UCSV')
fprintf('\n')
fprintf('Iteration 0000')

for irep = 1:nrep + nburn    % GIBBS iterations starts here
    % Print iterations
    if mod(irep,500) == 0
        fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%s%4d',8,8,8,8,8,8,8,8,8,8,8,8,8,8,'Iteration ',irep)
    end
    
    % Sample trend inflation (Bt) using Carter and Kohn (1994)    
    [Btdraw,log_lik] = carter_kohn_univ(y,Z,sigtS,sigtH,m,t,B_0_prmean,B_0_prvar);
    
    Btemp = Btdraw(:,2:t)' - Btdraw(:,1:t-1)';
    Bt2 = (Btemp).^2;
    bss = log(Bt2 + 1e-6);
               
    % Next draw statedraw (chi square approximation mixture component)
    % conditional on Sigtdraw       
    for jj=1:m   
        for i = 1:t-1
            for j = 1:numel(m_s)
                temp1= (1/sqrt(2*pi*u2_s(j)))*exp(-.5*(((bss(i,jj) - Htdraw(jj,i) - m_s(j) + 1.2704)^2)/u2_s(j)));
                prwH(j,1) = q_s(j,1)*temp1;       
            end
            prwH = prwH./sum(prwH);
            cprwH = cumsum(prwH);
            trand = rand(1,1);
            if trand < cprwH(1,1); imixH=1;
            elseif trand < cprwH(2,1), imixH=2;
            elseif trand < cprwH(3,1), imixH=3;
            elseif trand < cprwH(4,1), imixH=4;
            elseif trand < cprwH(5,1), imixH=5;
            elseif trand < cprwH(6,1), imixH=6;
            else imixH=7; 
            end
            statedrawH(i,jj)=imixH;   
        end
    end    
    
    vartH = u2_s(statedrawH);   % Variance
    bss1 = bss - m_s(statedrawH) + 1.2704;  % Mean 
    
    % Now sample the volatility of the state Bt (which I denote as Ht)
    % using CK
    [Htdraw] = carter_kohn(bss1',Zh,vartH,Qdraw,m,m,t-1,h_prmean,h_prvar);    
    % Final estimate of volatility
    sigtH = exp(Htdraw);
    Htdraw = Htdraw';
    
    % Qdraw is the variance of the state equation stochastic volatility
    sse_2 = (Htdraw(:,2:t-1) - Htdraw(:,1:t-2))*(Htdraw(:,2:t-1) - Htdraw(:,1:t-2))';
    Qinv = inv(sse_2 + Q_prmean);
    Qinvdraw = wish(Qinv,t+Q_prvar);
    Qdraw = inv(Qinvdraw);
    
    % Sample Sigma elements from Normal distribution
    y2 = [];
    for i = 1:t
        ytemps = y(:,i) - Z((i-1)*p+1:i*p,:)*Btdraw(:,i);
        y2 = [y2  (ytemps.^2)]; %#ok<AGROW>
    end
    yss = log(y2' + 1e-6);

    % Next draw statedraw (chi square approximation mixture component)
    % conditional on Sigtdraw       
    for i = 1:t
        for j = 1:numel(m_s)
            temp1= (1/sqrt(2*pi*u2_s(j)))*exp(-.5*(((yss(i,1) - Sigtdraw(1,i) - m_s(j) + 1.2704)^2)/u2_s(j)));
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
 
    %First draw volatilities conditional on sdraw
    vart = u2_s(statedrawS);
    yss1 = yss - m_s(statedrawS) + 1.2704;
    
    % Draw the measurement error volatility using CK
    [Sigtdraw] = carter_kohn(yss1',Zs,vart,Wdraw,1,1,t,sigma_prmean,sigma_prvar);
    sigtS = exp(Sigtdraw);
    
    Sigtdraw = Sigtdraw';
    % Draw the variance of the meaasurement error stochastic volatility
    sse_2 = (Sigtdraw(:,2:t) - Sigtdraw(:,1:t-1))*(Sigtdraw(:,2:t) - Sigtdraw(:,1:t-1))';
    Winv = inv(sse_2 + W_prmean);
    Winvdraw = wish(Winv,t+W_prvar);
    Wdraw = inv(Winvdraw);

    %----------------------------SAVE AFTER-BURN-IN DRAWS AND IMPULSE RESPONSES -----------------
    if irep > nburn;
        Bt_draws(:,irep-nburn) = Btdraw;
        Sigt_draws(:,irep-nburn) = sigtS;
    end % END saving after burn-in results 
end %END main Gibbs loop (for irep = 1:nrep+nburn)
fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8)
% fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8)
toc; % Stop timer and print total time
%=============================GIBBS SAMPLER ENDS HERE==================================
