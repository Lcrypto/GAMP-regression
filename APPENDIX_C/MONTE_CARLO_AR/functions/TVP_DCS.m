function [BETA_FINAL,SIGMA_FINAL] = TVP_DCS(y,x)

% **************************************************************************************************************************
% TVP_DCS.m  Code that estimates a univariate sparse, time-varying parameter (tvp) regression with constant variance. 
% The model is of the following form
%
%                   y[t]  =  x[t] beta[t] + e[t]
%                beta[t]  =  s[t] theta[t]
%               theta[t]  =  mu + 0.99I x (theta[t-1] - mu) + u[t]
%
% where e[t] ~ N(0, sigma), u[t] ~ N(0,(1/T)^2 x I), mu is the unconditional mean of the AR process for theta[t], s[t] is a 
% (px1) vector of 0/1 values (sparsity indicator), and we also assume the initial condition theta[0] ~ N(th0,V0).
% **************************************************************************************************************************
% Written by Dimitris Korobilis on 22/04/2017
% University of Essex
% **************************************************************************************************************************

[T,p] = size(x);

% Initial parameters   
lambda  = [0.001,0.01];      % prior probability of non-zero tap
p01     = [0.001,0.01];      % Pr{s_n(t) = 0 | s_n(t-1) = 1}
p10     = [0.001,0.01];      % Pr{s_n(t) = 1 | s_n(t-1) = 0}

eta     = [0,1,2];         % Mean of active coefficients
kappa   = 0.1;             % Circular variance of active coefficients
alpha   = [0.1,0.000001];  % Innovation rate of thetas (1 = total)     [dflt=0.10]
rho     = [1];            % Driving noise variance                    [dflt=(2 - alpha)*kappa/alpha]

imodel = 0;
nmodel = length(lambda)*length(p01)*length(p10)*length(alpha)*size(eta,2)*length(rho);
b_ols = [(x'*x)\(x'*y)];

PL = zeros(T,nmodel);

% Get index with all possible combinations of values
index = zeros(nmodel,6);
for i = 1:length(lambda)
    for j = 1:length(p01)
        for k = 1:length(p10)
            for l = 1:length(alpha)
                for m = 1:length(eta)
                    for n = 1:length(rho)
                        imodel = imodel + 1;
                        index(imodel,:) = [lambda(i),eta(m),alpha(l),rho(n),p01(j),p10(k)];
                    end
                end
            end
        end
    end
end

tic;
fprintf('Now you are running TVP_DCS')
fprintf('\n')
fprintf('\n')
fprintf('Iteration 00.00%% completed')
sigma = ones(T,1);
BETA_T = zeros(T,p,nmodel);
SIGMA_T = zeros(T,nmodel);
% Now estimate all possible models
for imodel=1:nmodel
    if mod(imodel,ceil(nmodel./10)) == 0
        fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%s%2.2f%s',8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,'Iteration ',100*(imodel/nmodel), '%% completed')
    end
    [MU_STATE,S_STATE] = DCS(y,x,sigma,index(imodel,1),index(imodel,2)*b_ols,kappa,index(imodel,3),index(imodel,4),index(imodel,5),index(imodel,6));
    if sum(sum(isnan(MU_STATE)))>0
        BETA_T(:,:,imodel) = zeros(T,p);
        SIGMA_T(:,imodel)  = ones(T,1);
        err = y - sum(x.*MU_STATE',2);
        PL(:,imodel) = normpdf(zeros(T,1),err,sqrt(sigma));
    else
        BETA_T(:,:,imodel) = MU_STATE';
        SIGMA_T(:,imodel)  = S_STATE;
        err = y - sum(x.*MU_STATE',2);
        PL(:,imodel) = 0*normpdf(zeros(T,1),err,sqrt(sigma)) + 1e-20;
    end

end

probs = PL./repmat(sum(PL,2),1,nmodel);
BETA_FINAL = 0*MU_STATE';
SIGMA_FINAL = 0*S_STATE;
for i = 1:nmodel
    BETA_FINAL = BETA_FINAL + repmat(probs(:,i),1,p).*BETA_T(:,:,i);
    SIGMA_FINAL = SIGMA_FINAL + probs(:,i).*SIGMA_T(:,i);
end

fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8)
toc;