% KP_PPT.m model which does structural breaks using KP(2007) prior in
% conditional mean and variance, and PPT(2006) prior for the break process.
% -------------------------------------------------------------------------
% The model is:
%
%         y[t] = Z[t] beta[s[t]] + exp( sigma[s[t]] / 2) e[t]
%         beta[m] = beta[m-1] + u[t]
%         sigma[m] = sigma[m-1] + v[t]
% 
%  ...(add later)
% -------------------------------------------------------------------------

clear all;
clc;
tic;

% ==========================| LOAD DATA
% Load the different series
load GDP.dat;
load HOUSE.dat;
load IP.dat;
load M2.dat
load PCE.dat;
load TBILL.dat;
%---------------
% Choose one of the series above
Yraw=GDP(1:end-8,:);

%===========================| USER INPUT |=================================
% Gibbs related preliminaries
nrep = 1000;           % Gibbs draws to keep
nburn = 1000;          % Burn in phase
ntot = nrep + nburn;  % Total draws
it_print = 50;        % Print every "it_print"-th iteration

% Set stationarity transformation for dependent variable y
tcode = 1;            % 1: levels
                      % 2: first differences
                      % 3: second differences
                      % 4: logs
                      % 5: first log differences
                      % 6: second log differences

% Model specification
constant = 1;         % 0: no intercept; 1: intercept
lags = 1;             % Number of AR lags

M = 3;                % Number of MAXIMUM break points
                      % Implies that there is a MAXIMUM M+1 periods

forecasting = 0;      % 0: no forecasting; 1: forecasting
nfore = 8;            % Number of predictions steps ahead
nforerep = 10;       % Nubmer of times to forecast at each MCMC iteration

%===========================| PRELIMINARIES |==============================
% ----Transform to stationarity, if applicable
% First evaluate user input for tcode
if tcode<1 && tcode>6
    error('Wrong choice of transformation code (tcode)');
end
ytempraw = transx(Yraw,tcode);
% For differenced data, we loose one observation
if tcode ~= 1 && tcode ~= 4
    ytempraw = ytempraw(2:end,:);
end

% Get time series observations
traw = size(ytempraw,1);

% Create RHS variables, x
ylag = mlag2(ytempraw,lags);  % This routine creates the lagged dep variables   
if constant == 1
    X = [ones(traw-lags,1) ylag(lags+1:end,:)];
    k = lags + 1;
elseif constant == 0
    X = ylag(lags+1:end,:);
    k = lags;
else
    error('Wrong specification of intercept')
end

% Correct # of observations of y after taking lags
y = ytempraw(lags+1:end,:);
t = size(y,1);

% Get OLS quantities from AR(lags) model
res=ols(y,X);
disp('');
disp('Ols results');
disp('');
prt(res);
beta_hat = y'/X';
sigma_hat = (y-X*beta_hat')'*(y-X*beta_hat')./(t-k-1);

%===========================| PRIORS
% Transition matrix probabilities
% pii distr as a Beta(a,b)
a=1; b=1;                       %1. Non informative priors
%a=0.5*sqrt(t); b=0.5*sqrt(t);   %2. Least informative Beta prior
%a=0.5; b=0.5;                   %3. Reference prior
%a=0.1; b=10;                      %4. Informative prior (few breaks);

% For sigma state equation innovation precision, G(_eta1,1/_eta2)
eta1 = .01;
eta2 = .01;
eta = (eta2/eta1);  % Prior mean of eta (as initial value)

% Hyperparameters for beta state equation innovation, W(hyp_sdf,scalmat) 
hyp_sdf = (k+1);
scalmat = 1*eye(k);
V = 0.1*eye(k,k); % Initialize V, the state variance of betas

%===========================| STORAGE SPACE FOR POSTERIORS & STUFF
% Transition probability matrix
P(:,:)=0.5.*eye(M+1,M+1);
for i=1:M       
    P(i,i+1)=1-P(i,i);
end

% Initialize betas, sigma2
betas = 0.2*ones(1,k*(M+1));  %repmat(beta_hat,1,M+1);
sigma2 = 0.1*ones(1,M+1);   %repmat(sigma_hat,1,M+1); 

nu1 = 0.01*ones(t,1);
nu2 = nu1;
Bdraw = ones(t,k);

betasdraws = zeros(nrep,k*(M+1));
sigma2draws = zeros(nrep,M+1);
Vdraws = zeros(nrep,k,k);
etadraws = zeros(nrep,1);
sdraws = zeros(nrep,t);
Pdraws = zeros(nrep,M+1,M+1);
y_foredraws = zeros(nrep*nforerep,nfore);

%==========================================================================
%==========================| GIBBS SAMPLER |===============================
for irep = 1:ntot    % GIBBS iterations start here
    % Print iterations
    if mod(irep,it_print) == 0
        disp(irep);toc;
    end

    %-----| STEP 1: DRAW STATES, S[T]
    [pst,s,llikf] = drawstatesun(y,X,M,betas,sigma2,P);   

    %-----| STEP 2: DRAW PARAMETERS, BETAS, SIGMA2
    % First get durations in each regime:
    dur = zeros(M+1,1);
    dur(1,1)=1;
    for i=1:M+1
        dur(i,1) = sum(s==i);
    end
    
    M_new = M;%s(t) - 1;
    
    % restrict very small regimes
    s1=s;
    prev_dur = 0;
    for i = 1:M+1
        prev_dur = sum(dur(1:i-1));
        if dur(i,1)==1
            if prev_dur+dur(i)<t
                s1(prev_dur+1:prev_dur+dur(i,1)) = s(prev_dur+dur(i,1)+1);
                s1(prev_dur+1:end) = s1(prev_dur+1:end)-1;
            else %we are in the last regime
                s1(prev_dur+1:prev_dur+dur(i,1)) = s(prev_dur);
            end
        end
    end
    s=s1;
    
    % 2.i) Sample regression coefficients, betas---------------------------
    t_m = zeros(t,k);  %storage matrix for state mean
    t_v = zeros(k,t*k);  %storage matrix for state variance
    c_m = zeros(k,1);   %initial state mean is 0
    t_m(1,:) = c_m';
    v_m = 100*eye(k);   %initial state variance is 100
    t_v(:,1:k) = v_m;
    for i = 1:t
        h = X(i,:);     
        pst_cf = c_m;   %Predicted state mean   
        if i>1 && s(i-1)~=s(i);
            %if there is switch between t and t+1:       
            pst_vr = v_m + V; %Predicted state covariance       
        else
            %no switch - do not update predicted covariance with sigma
            pst_vr = v_m; %Predicted state covariance
        end
        % Kalman filter iterations
        f_cast = y(i,:) - h*pst_cf;  %cond. forecast error
        ss = h*pst_vr*h' + nu1(i,1);  
        K_GN = pst_vr*h'*inv(ss); %Kalman gain
        c_m = pst_cf+K_GN*f_cast;  %updated state mean
        v_m = (eye(k)-K_GN*h)*pst_vr;   %updated state variance
        t_m(i,:) = c_m';     %store state mean
        t_v(:,(i-1)*k+1:i*k) = v_m;  %store state variance  
    end
    
    tt_m = t_m;
    tt_v = t_v;
    % Draw the last (at time T) theta, using mean and variance from forward iterations above
    Bdraw(s(t),:) = tt_m(t,:) + (chol(tt_v(:,(t-1)*k+1:t*k))'*randn(k,1))'; %#ok<*AGROW>

    if s(t,:)<t;       
        for i = s(t)+1:t
            Bdraw(i,:) = Bdraw(i-1,:) + (chol(V)'*randn(k,1))';
            vv((i-2)*k+1:(i-1)*k,:) = (Bdraw(i,:)-Bdraw(i-1,:))'*(Bdraw(i,:)-Bdraw(i-1,:));
            nu2(i,:) = exp(log(nu2(i-1,:)) + sqrt(eta)*randn);
            vvv(i-1,:) = (log(nu2(i,:)) - log(nu2(i-1,:)))^2;
                
        end
    end
    
    %start backward iterations
    for i = t:-1:2       
        if s(i) ~= s(i-1)
            f_cast0m = Bdraw(s(i),:)-tt_m(i-1,:);       
            ss0m = tt_v(:,(i-2)*k+1:(i-1)*k) + V;       
            K_GN0m = tt_v(:,(i-2)*k+1:(i-1)*k)*inv(ss0m);
            tt_m(i-1,:) = tt_m(i-1,:) + f_cast0m*K_GN0m;               
            tt_v(:,(i-2)*k+1:(i-1)*k) = (eye(k)-K_GN0m)*tt_v(:,(i-2)*k+1:(i-1)*k);
            Bdraw(s(i-1),:) = tt_m(i-1,:) + (chol(tt_v(:,(i-2)*k+1:(i-1)*k))'*randn(k,1))';
            vv((s(i-1)-1)*k+1:s(i-1)*k,:) = (Bdraw(s(i),:)-Bdraw(s(i-1),:))'*(Bdraw(s(i),:)-Bdraw(s(i-1),:));
        end
    end
    
    % Btdraw is the Tx1 vector of parameters that will be used to condition
    % on when estimating the volatilities Sigtdraw
    Btdraw = [];
    for i = 1:M+1      
        Btdraw = [Btdraw; repmat(Bdraw(i,:),dur(i,1),1)]; %#ok<AGROW>
    end
    
    y2 = [];
    for i = 1:t
        ytemps = y(i,:) - X(i,:)*Btdraw(i,:)';   
        y2 = [y2  (ytemps.^2)]; %#ok<AGROW>
    end
    
    sigma2(1,1) = sum(y2(1:dur(1,1)))/(dur(1,1));
    prev_dur = 0; %#ok<NASGU>
    for i = 2:M+1
        prev_dur = sum(dur(1:i-1));   
        if dur(i)~=0
            sigma2(1,i) = sum( y2(:,prev_dur+1:prev_dur+dur(i)) )/(dur(i));
        else
            sigma2(1,i) = sigma2(1,i-1);
        end
        sigma2(sigma2(1,i)<0.000001,i) = sigma2(1,i-1);
    end
    
    nu2(1,:) = sigma2(1,1);
    vvv(1,:) = (log(sigma2(:,1)))^2;
    for i = M:-1:1
        nu2(i+1,:) = sigma2(:,i+1);
        vvv(i,:) = (log(sigma2(:,i+1)) - log(sigma2(:,i)))^2;
    end
    
    % 2.ii) Sample the state covariances V,eta-----------------------------
    schol2 = scalmat; %prior scale matrix for sigma
    sum3 = eta2;      %prior scale matrix for eta
    for i = 1:t-1
        schol2 = schol2 + vv((i-1)*k+1:i*k,:); %sse for theta state equation
        sum3 = sum3 + 0.5*vvv(i,:);             %sse for nu state equation
    end
    Vinv = inv(schol2);
    Vinvdraw = wish(Vinv,t+hyp_sdf);          %draw V^-1
    V = inv(Vinvdraw);
    eta = inv( gamrnd(eta1+.5*t,1/sum3) );      %draw eta
    if eta > 5
        eta=1;
    end 
    
    betas = zeros(1,(M+1)*k);
    betas(1,1:(M+1)*k)=reshape(Bdraw(1:M+1,:)',1,(M+1)*k);
    
    sigma2 = zeros(1,M+1);
    sigma2(1,1:M_new+1) = nu2(1:M+1,:)';
    
    % nu is the Tx1 vector of volatilities that will be used to condition
    % on when estimating the Bdraws
    nu1 = [];
    for i = 1:M+1
        nu1 = [nu1; repmat(sigma2(1,i),dur(i,1),1)]; %#ok<AGROW>
    end
    
    %-----| STEP 3: DRAW TRANSITION PROBABILITIES, P
    P = drawPun(y,X,M,s,a,b);   
    
    %================Save post burn-in draws
    if irep>nburn
        %save parameter draws
        betasdraws(irep-nburn,:) = betas;
        sigma2draws(irep-nburn,:) = sigma2;
        etadraws(irep-nburn,:) = eta;
        Vdraws(irep-nburn,:,:) = V;        
        sdraws(irep-nburn,:) = s';
        Pdraws(irep-nburn,:,:) = P;
        
        %==========Forecasting
        if forecasting == 1
            y_fore = zeros(nforerep,nfore);
            for forep = 1:nforerep
                s_fore = [s ; zeros(nfore,1)];
                betas_fore = [betas(1,(s(t)-1)*k+1:s(t)*k) ; zeros(nfore,k)];
                sigma2_fore = [sigma2(1,s(t)); zeros(nfore,1)];
                for ii = 1:nfore
                    nii = sum(s_fore==s_fore(t+ii-1));
                    p_fore = betarnd(a+nii,b+1);    % Draw Pkk
                    % Draw S(T+h)
                    if p_fore > rand
                        s_fore(t+ii) = s_fore(t+ii-1);       
                    else
                        s_fore(t+ii) = s_fore(t+ii-1) + 1;
                    end
                    
                    % Update the parameters
                    if s_fore(t+ii) == s_fore(t+ii-1)
                        betas_fore(ii+1,:) = betas_fore(ii,:);
                        sigma2_fore(ii+1,:) = sigma2_fore(ii,:);
                    else
                        betas_fore(ii+1,:) = betas_fore(ii,:) + randn(1,k)*chol(V);
                        % We do not want to forecast with explosive
                        % parameter draws. Here I check only the AR(1)
                        % coefficient and if it is explosive I set it back
                        % to the previous regime (instead of drawing again
                        % until stationary which can be inefficient computationally)
                        if abs(betas_fore(ii+1,1*constant + 1)) >= 1.000                        
                           betas_fore(ii+1,1*constant + 1) = betas_fore(ii,1*constant + 1);
                        end
                        log_sigma_tplusone = log(sigma2_fore(ii,:)) + sqrt(eta)*randn;
                        sigma2_fore(ii+1,:) = exp(log_sigma_tplusone);
                    end
                end
                X_fore = [1 y(t,:) X(t,2:lags)];
                y_hat = X_fore*betas_fore(2,:)' + sqrt(sigma2_fore(2))*randn;
                y_fore(forep,1) = y_hat;
                for ii = 1:nfore-1  % Predict T+2, T+3 until T+h                             
                    if ii < lags 
                        X_fore = [1 y_hat X_fore(:,2:lags)];                             
                        % This gives the forecast T+i for i=1,..,p                             
                        y_hat = X_fore*betas_fore(ii+2,:)' + sqrt(sigma2_fore(ii+2))*randn;     
                        y_fore(forep,ii+1) = y_hat;
                    elseif  ii >= lags
                        X_fore = [];       
                        for indlag = 1:lags                  
                            X_fore =  [X_fore y_fore(forep,ii - indlag + 1)];        %#ok<AGROW>
                        end
                        X_fore = [1 X_fore]; %#ok<AGROW>
                        y_hat = X_fore*betas_fore(ii+2,:)' + sqrt(sigma2_fore(ii+2))*randn;
                        y_fore(forep,ii+1) = y_hat;
                    end
                end %  the last value of 'Y_temp' is the prediction T+h
                
            end  % End forecasting nforerep times
            y_foredraws(((irep-nburn)-1)*nforerep+1:(irep-nburn)*nforerep,:) = y_fore;
        end  % End forecasting 
    end
    
    
end


