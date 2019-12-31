
function [storage_Bt,storage_SIGMAt] = TVP_LASSO(y,z,LASSO,nsave,nburn)
%==========================================================================
% Code to replicate results for Belmonte, Koop & Korobilis (2011)
%--------------------------------------------------------------------------
% The final state-space form of the model is:
%
%            y[t+h]  = b x[t] + omega b*[t] x[t] + e[t+h]
%            b*[t]   = b*[t-1] + u[t]
% 
% where e[t]~N(0,sigma2), u[t]~N(0,1) and b*[t] = b[t]/omega.
%--------------------------------------------------------------------------
% OUTPUT:
%      Bt      Draws of TVP regression coefficients
%      St      Draws of stochastic variances
% INPUT:
%       y      Dependent variable
%       z      Predictor variables
%   LASSO      1: LASSO everywhere; 2: LASSO only on the initial condition
%              3: LASSO only on the state variance omega; 4: No LASSO
%--------------------------------------------------------------------------
% Written by Dimitris Korobilis, 2011
% CORE, Universite Catholiue de Louvain
%==========================================================================


[T,k] = size(z);
ntot = nsave + nburn;

%========= PRIORS
% Prior on the LASSO regularization parameters, lambda2 (for the initial condition):
a1_l = 3; a2_l = 1;
% Prior on the LASSO regularization parameters, kappa2 (for the variance of the state):
a1_k = 3; a2_k = 1;

% Prior on stochastic volatility:   
% Initial condition, log(sigma_0) ~ N(log(sigma_OLS),I_n)   
sigma_prmean = 0;
sigma_prvar = 10;
    
% Volaitlity variance, W ~ IW(k2_W*(1+dimension(W))*I_n,(1+dimension(W)))
W_prmean = 1;   
W_prvar = 12;

% Parameters of the 7 component mixture approximation to a log(chi^2)
% density:
q_s = [0.00730; 0.10556; 0.00002; 0.04395; 0.34001; 0.24566; 0.25750];
m_s = [-10.12999; -3.97281; -8.56686; 2.77786; 0.61942; 1.79518; -1.08819];
u2_s = [5.79596; 2.61369; 5.17950; 0.16735; 0.64009; 0.34023; 1.26261];
offset = 1e-6;

%========= INITIALIZE MATRICES:
lambda2 = 1;
tau2 = 2*ones(k,1);
kappa2 = 1;
ksi2 = 2*ones(k,1);
omega = 0.2*ones(k,1);
smoothed_states = 0.01*ones(T,k);
XB = ( ( z .* smoothed_states) * omega ); 

% Prior variance on init.cond. (b0) and state variance (omega) when NOT using the LASSO prior
Vb0_inv = inv(9*eye(k));
Vomega_inv = inv(4*eye(k));

% Initialize matrices for stochastic volatility equation
consW = 0.01;
Wdraw = consW;  % This is the variance of the stochastic volatility
Sigtdraw = 0.1*ones(T,1);
sigtS = 0.2*ones(T,1);
statedrawS = 5*ones(T,1);
Zs = ones(T,1);
prwS = zeros(numel(q_s),1);

% Storage matrices for posteriors and stuff
storage_Bt = zeros(T,k,nsave);
storage_SIGMAt = zeros(T,nsave);

%============================= START SAMPLING ===================================
%================================================================================
tic;
fprintf('Now you are running TVP_LASSO')
fprintf('\n')
fprintf('Iteration 0000')

for irep = 1:ntot
        % Print iterations   
        if mod(irep,500) == 0
            fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%s%4d',8,8,8,8,8,8,8,8,8,8,8,8,8,8,'Iteration ',irep)
        end
        
        temp = 1./sigtS;
        inv_sigma2 = diag(temp); %This is the inverse of the SV: ( exp(h[t]/2) )^2
        
        %% STEP I: sampling initial condition b0
        y_tilde = y - XB;        
        Dinv =  ((z'*inv_sigma2*z) + Vb0_inv)\eye(k); 
        Sigma_b0 = Dinv ;
        mean_b0 = Dinv*(z'*inv_sigma2*y_tilde); 
        b0 = mean_b0 + chol(Sigma_b0)'*randn(k,1);

        if LASSO==1 || LASSO==2
            %% STEP II: Sample tau2 and update Cinv
            mean_tau2 = sqrt(lambda2./(b0.^2));
            var_tau2 = lambda2.*ones(k,1);
            inv_tau2 = inv_gaussian_generator(mean_tau2,var_tau2,k) +1e-8;
            Vb0_inv = diag(inv_tau2) ;

            %% STEP III: Sample lambda2 (lasso parameter)
            shape_t = k + a1_l;
            tau2 = 1./inv_tau2; 
            scale_t = 0.5*sum(tau2) + a2_l;
            lambda2 = gamrnd(shape_t,1/scale_t, 1,1);
        end
        
        %% STEP IV: sample omegas
        centered = y - (z*b0);
        Z = (z.*smoothed_states);
        Ainv = ((Z'*inv_sigma2*Z) + Vomega_inv)\eye(k); 
        Sigma_omega = Ainv ; 
        mean_omega = Ainv *(Z'*inv_sigma2*centered);
        omega = mean_omega + chol(Sigma_omega)'*randn(k, 1);

        if LASSO==1 || LASSO==3
            %% STEP V: sample epsilon2
            mean_ksi2 = sqrt(kappa2./(omega.^2));       
            var_ksi2 = kappa2.*ones(k,1);
            inv_ksi2 = inv_gaussian_generator(mean_ksi2,var_ksi2,k) +1e-8;
            Vomega_inv = diag(inv_ksi2);
 
            %% STEP VI: Sample kappa2
            shape_k = k + a1_k;
            ksi2 = 1./inv_ksi2; 
            scale_k = 0.5*sum(ksi2) + a2_k;
            kappa2 = gamrnd(shape_k,1/scale_k, 1,1);
        end
        
        %% STEP VII: Now update tvp with updated omega and sigma2 parameters
        Z = (z.*repmat(omega',T,1)) ;
        Intercept = (z*b0);
        [smoothed_states,~] = carter_kohn((y-Intercept)',Z,sigtS,eye(k),k,1,T,zeros(k,1),zeros(k,k));
        
        s1 = rand(1,1);
        if s1>0.5 % we implement a switch in sign with probability 0.5 in omegas and regressors
            omega = -omega;
            smoothed_states = -smoothed_states;
        end
        
        % update regressors and residuals for next iteration.
        XB = ((z .* smoothed_states)*omega); 
        e = (y - XB - Intercept);
        
        %% STEP VII: Now update stochastic volatility in the measurement
        % equation (volatility of inflation)
        yss = log(e.^2 + 1e-6);
        
        % Next draw statedraw (chi square approximation mixture component)   
        % conditional on Sigtdraw       
        for i = 1:T
            for j = 1:numel(m_s)
                temp1= (1/sqrt(2*pi*u2_s(j)))*exp(-.5*(((yss(i,1) - Sigtdraw(i,1) - m_s(j) + 1.2704)^2)/u2_s(j)));
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
        [Sigtdraw,~] = carter_kohn(yss1',Zs,vart,Wdraw,1,1,T,sigma_prmean,sigma_prvar);
        
        sigtS = exp(0.5*Sigtdraw);
        
        Sigttemp = Sigtdraw(2:T,:) - Sigtdraw(1:T-1,:);   
        sse_2 = sum(Sigttemp.^2);
                
        % Draw the variance of the meaasurement error stochastic volatility
        W = 0.5*sse_2 + W_prmean;
        Winvdraw = gamrnd(W_prvar+.5*T,1/W);
        Wdraw = inv(Winvdraw);
        
        %================Save post burn-in draws
        if irep>nburn
            storage_Bt(:,:,irep-nburn) = repmat(b0',T,1) + smoothed_states;
            storage_SIGMAt(:,irep-nburn) = sigtS;
        end
end
fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8)
% fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8)
toc; % Stop timer and print total time

