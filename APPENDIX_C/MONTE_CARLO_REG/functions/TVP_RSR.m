function [y_t_DMA,var_DMA,log_PL_DMA]  = TVP_RSR(y_t,X,X_fore,plag)
%==========================================================================
% TVPBCR.m: Time-varying parameter Bayesian compressed regression for
%           forecasting inflation, as in Korobilis and Pettenuzzo (2016)
%==========================================================================
% TVP-AR with forgetting and recurssive moment estimation of the
% measurement variance.
%
%        y[t] = theta[t] x Q x z[t] + e1,        e1 ~ N(0,V_t)
%    theta[t] = theta[t-1]      + e2,        e2 ~ N(0,S_t)  
%
% where z[t] is the matrix of predictor variables, theta[t] is the
% time-varying regression coefficient, and Q is a random projection matrix.
%==========================================================================
% Written on 20/08/2015
% Dimitris Korobilis,
% University of Glasgow
%==========================================================================


% =============================| MODEL SPECIFICATION |========================= 
% Forgetting factors
lambda = 0.99;            % For the time-varying parameters theta
alpha = 0.98;             % For the model switching
kappa = 0.94;             % For the error covariance matrix

                          
% Initial values on time-varying parameters
% theta[0] ~ N(PRM,PRV x I)
prior_theta = 2;          % 1: Diffuse N(0,4)
                          % 2: Data-based prior
% =============================| end model specification |=========================


%=================================| PRELIMINARIES |================================
Z_t = X(:,1:1+plag);
z_t = X(:,2+plag:end);

Z_t_fore = X_fore(:,1:1+plag);
z_t_fore = X_fore(:,2+plag:end);

[T,m] = size(X);

%============| DEFINE RANDOM PROJECTION:
% We need to define basics of random projections. Number of RHS variables is m (including intercept and "plag" lags of the dependent)
m_max = max(m/2,5*round(log(m)));      % This is the max size. Note that m_max<<m
m_min = 1*round(log(m));             % This is the min model size
RP_type = 6;                % Type of random projection (see function getRP.m)
n_psi = 10;                 % Times to simulate random projection for a given model size

index_models = repmat((m_min:m_max)',n_psi,1);
K = length(index_models);

% Get size of projected variables
m_bcr = size(z_t,2);

% Initial value of measurement error covariance matrix V_t
V_0 = var(y_t(1:20,:))./4;
prob_0_prmean = 1./K;

% Define forgetting factor(s):                 
inv_lambda = 1./lambda;
expert_weight = 1;

% Initialize matrices in PC memory
theta_0_prmean = cell(K,1);
theta_0_prvar = cell(K,1);
theta_pred = cell(K,1);
R_t = cell(K,1);
prob_pred = zeros(T,K);
y_t_pred = cell(K,1);
y_t_pred_h = cell(K,1);
e_t = cell(K,1);
A_t = cell(K,1);
V_t = cell(K,1);
theta_update = cell(K,1);
S_t = cell(K,1);
variance = cell(K,1);
w_t = cell(K,1);
log_PL = zeros(K,1);
prob_update = zeros(T,K);
y_t_DMA = zeros(T,1);
var_DMA = zeros(T,1);
y_t_BEST = zeros(T,1);
var_BEST = zeros(T,1);
log_PL_DMA = zeros(T,1);
log_PL_BEST = zeros(T,1);
xRx = cell(K,1);
offset = 1e-20; % This offset constant is used in some cases for numerical stability   
% =============================| end of preliminaries |=========================
tic;
fprintf('Now you are running TVP_RSR')
fprintf('\n')
fprintf('\n')
fprintf('Iteration 00.00%% completed')
% =============================Start now the Kalman filter loop   
for irep = 1:T % for 1 to T time periods
    if mod(irep,ceil(T./10)) == 0
        fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%s%2.2f%s',8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,'Iteration ',100*(irep/T), '%% completed')
    end
    
    % Here get the sum of all K model probabilities, quantity you
    % are going to use to update the individual model probabilities
    if irep>1           
        % Exponential Forgetting
        sum_prob_a = sum((prob_update(irep-1,:).^alpha).*(expert_weight^(1-alpha)),2);  % this is the sum of the K model probabilities (all in the power of the forgetting factor 'a')
    end

    % reset log_PL, A_t and R_t, to zero at each iteration to save memory
    log_PL = zeros(K,1);    
    A_t = cell(K,1);
    R_t = cell(K,1);
    for k = 1:K % for 1 to K projected models
        % Generate PHI
        [PHI] = genRP(RP_type,index_models(k),m_bcr);
        
        x_t = cell(1,1);  x_t_fore = cell(1,1);
        x_t{1,1} = [Z_t z_t*PHI']; %#ok<*USENS>
        x_t_fore{1,1} =  [Z_t_fore z_t_fore*PHI'];
        % -----------------------Predict
        if irep==1
            %%%%%%%%%%%%%%%%%%%
            % Prior is defined here for BCR
            m_tot = size(x_t{1,1},2);
            if prior_theta == 1 % diffuse
                theta_0_prmean{ll,1} = 0*ones(m_tot,1);
                theta_0_prvar{ll,1} =  4*eye(m_tot);
            elseif prior_theta == 2 % data-based 1
                varx = var(x_t{1,1});
                varx(varx==0)=0.01;
                vary = var(y_t(1:20,:));                
                theta_0_prmean{k,1} = 0*ones(m_tot,1);
                theta_0_prvar{k,1} =  2*diag((vary./varx)');
            end
            %%%%%%%%%%%%%%%%%%%
            theta_pred{k,1} = theta_0_prmean{k,1};  % predict theta[t], this is Eq. (5)
            R_t{k,1} = inv_lambda*theta_0_prvar{k,1};   % predict R[t], this is Eq. (6)
            temp1 = ((prob_0_prmean).^alpha);  
            prob_pred(irep,k) = temp1./(K*temp1);     % predict model probability, this is Eq. (15)
        else
            theta_pred{k,1} = theta_update{k,1};    % predict theta[t], this is Eq. (5)
            R_t{k,1} = inv_lambda.*S_t{k,1};   % predict R[t], this is Eq. (6)
            % Exponential Forgetting           
            prob_pred(irep,k) = ((prob_update(irep-1,k).^alpha)*(expert_weight^(1-alpha)) + offset)...
                ./(sum_prob_a + offset);   % predict model probability, this is Eq. (15)            
        end

        % Now implememnt individual-model predictions of the variable of interest
        y_t_pred{k,1}(irep,:) = x_t{1,1}(irep,:)*theta_pred{k,1};   %one step ahead prediction
        
        % Now do h_fore-step ahead prediction
        y_t_pred_h{k,1}(irep,:) = x_t_fore{1,1}(1,:)*theta_pred{k,1}; % predict t+h given t
        
        % -------------------------Update
        e_t{k,1}(:,irep) = y_t(irep,:) - y_t_pred{k,1}(irep,:); % one-step ahead prediction error
        
        % We will need some products of matrices several times, which is better to define them
        % once here for computational efficiency
        R_mat = R_t{k,1};
        xRx2 = x_t{1,1}(irep,:)*R_mat*x_t{1,1}(irep,:)';
        
        % Update V_t - measurement error covariance matrix using rolling
        % moments estimator, see top of page 12
        if irep==1
            V_t{k,1}(:,irep) = V_0;
        else
            A_t{k,1} = (e_t{k,1}(:,irep-1)).^2;
            V_t{k,1}(:,irep) = kappa*V_t{k,1}(:,irep-1) + (1-kappa)*A_t{k,1};
        end
        
        % Update theta[t] (regression coefficient) and its covariance
        % matrix S[t]
        Rx = R_mat*x_t{1,1}(irep,:)';
        KV = V_t{k,1}(:,irep) + xRx2;
        KG = Rx/KV;
        theta_update{k,1} = theta_pred{k,1} + KG*e_t{k,1}(:,irep);
        S_t{k,1} = R_mat - KG*(x_t{1,1}(irep,:)*R_mat); %#ok<*MINV>
        
        % Update model probability. Feed in the forecast mean and forecast
        % variance and evaluate at the future inflation value a Normal density.
        % This density is called the predictive likelihood (or posterior
        % marginal likelihood). Call this f_l, and use that to update model
        % weight/probability called w_t
        variance{k,1}(irep,:) = V_t{k,1}(:,irep) + xRx2;   % This is the forecast variance of each model
        if variance{k,1}(irep,:)<=0  % Sometimes, the x[t]*R[t]*x[t]' quantity might be negative
            variance{k,1}(irep,:) = abs(variance{k,1}(irep,:));
        end
        mean_f = x_t{1,1}(irep,:)*theta_pred{k,1};  % This is the forecast mean
        f_l = (1/sqrt(2*pi*variance{k,1}(irep,:)))*exp(-.5*(((y_t(irep,:) - mean_f)^2)/variance{k,1}(irep,:))); %normpdf(y_t(irep,:),mean_f,sqrt(variance));
        w_t{k,1}(:,irep) = prob_pred(irep,k)*f_l;
        
        % Calculate log predictive likelihood for each model
        log_PL(k,1) = log(f_l + offset);
    end % end cycling through all possible K models
    
    % First calculate the denominator of Equation (16) (the sum of the w's)
    sum_w_t = 0;
    for k_2=1:K %#ok<*BDSCI>
        sum_w_t = sum_w_t + w_t{k_2,1}(:,irep);
    end
    
    % Then calculate the updated model probabilities
    for k_3 = 1:K
        prob_update(irep,k_3) = (w_t{k_3,1}(:,irep) + offset)./(sum_w_t + offset);  % this is Equation (16)
    end
    
    if irep == T
        % Now we have the predictions for each model & the associated model probabilities: Do DMA forecasting
        for k_4 = 1:K
            model_i_weight = prob_pred(irep,k_4);
            % The next temp_XXX calculate individual model quantities, weighted by their model probabilities. Then take the sum of all these.
            temp_pred = y_t_pred_h{k_4,1}(irep,:)*model_i_weight;       
            temp_var = variance{k_4,1}(irep,:)*model_i_weight;
            temp_logPL = log_PL(k_4,1)*model_i_weight;
            y_t_DMA(irep,:) = y_t_DMA(irep,:) + temp_pred;  % This is the mean DMA forecast
            var_DMA(irep,:) = var_DMA(irep,:) + temp_var;   % This is the variance of the DMA forecast
            log_PL_DMA(irep,:) = log_PL_DMA(irep,:) + temp_logPL;  % This is the DMA Predictive Likelihood   
        end
    end
    
    % Get log_PL_BEST here (cannot get it after the main loop is finished, like with y_t_BEST)
    [temp_max_prob, temp_best_model] = max(prob_update(irep,:));
    log_PL_BEST(irep,:) = log_PL(temp_best_model,:);
end
%***********************************************************
% Find now the best models 
max_prob = zeros(T,1);
best_model = zeros(T,1);
for ii=1:T
    [max_prob(ii,1), best_model(ii,1)]=max(prob_pred(ii,:));
end

% ...and make the prediction based on the best model
for ii=1:T
    y_t_BEST(ii,1) = y_t_pred_h{best_model(ii,1),1}(ii,:);
    var_BEST(ii,1) = variance{best_model(ii,1),1}(ii,:);
end

fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8)
toc;


