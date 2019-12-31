function [y_foredraws] = KP_PPT(y,X,X_fore,constant,nrep,nburn)

% KP_PPT.m model which does structural breaks using KP(2007) prior in
% conditional mean and variance, and PPT(2006) prior for the break process.

M = 2;

[t,k] = size(X);

% Get OLS quantities
beta_hat = y'/X';
sigma_hat = (y-X*beta_hat')'*(y-X*beta_hat')./(t-k-1);

%===========================| PRIORS
% Transition matrix probabilities
% pii distr as a Beta(a,b)
a=1; b=1;

% For sigma state equation innovation precision, G(_eta1,1/_eta2)
eta1 = .01;
eta2 = .01;
eta = (eta2/eta1);  % Prior mean of eta (as initial value)

% Hyperparameters for beta state equation innovation, W(hyp_sdf,scalmat) 
hyp_sdf = (k+2);
scalmat = 1*eye(k);
V = 0.01*eye(k,k); % Initialize V, the state variance of betas

%===========================| STORAGE SPACE FOR POSTERIORS & STUFF
% Transition probability matrix
P(:,:)=0.5.*eye(M+1,M+1);
for i=1:M       
    P(i,i+1)=1-P(i,i);
end

% Initialize betas, sigma2
betas = 0.3*ones(1,k*(M+1));  %repmat(beta_hat,1,M+1);
sigma2 = 0.1*ones(1,M+1);   %repmat(sigma_hat,1,M+1); 
sigma_old = sigma2;

nu1 = 0.01*ones(t,1);
nu = nu1;
Bdraw = ones(t,k);

betasdraws = zeros(nrep,k*(M+1));
sigma2draws = zeros(nrep,M+1);
Vdraws = zeros(nrep,k,k);
etadraws = zeros(nrep,1);
sdraws = zeros(nrep,t);
Pdraws = zeros(nrep,M+1,M+1);
y_foredraws = zeros(nrep,1);
%==========================================================================
%==========================| GIBBS SAMPLER |===============================
tic;
fprintf('Now you are running KP_PPT')
fprintf('\n')
fprintf('Iteration 0000')
for irep = 1:nrep+nburn    % GIBBS iterations start here
    % Print iterations
    if mod(irep,500) == 0
        fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%s%4d',8,8,8,8,8,8,8,8,8,8,8,8,8,8,'Iteration ',irep)
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

    % restrict regimes with duration only 1 period
    s1=s;
    prev_dur = 0;
    for i = 1:M+1
        prev_dur = sum(dur(1:i-1));
        if dur(i,1)<=3
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
    [Bdraw,vv,nu,vvv] = drawpars(y,X,nu1,nu,V,eta,t,k,M,s,dur);
    
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
    if eta > 8
        eta = .1;
    end
    
    betas = zeros(1,(M+1)*k);
    betas(1,1:(M+1)*k)=reshape(Bdraw(1:M+1,:)',1,(M+1)*k);
        
    sigma2 = zeros(1,M+1);
    sigma2(1,1:M+1) = nu(1:M+1,:)';
    sigma2(sigma2<0.0001)=sigma_old(sigma2<0.0001);
    sigma_old = sigma2;
    
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
        betasdraws(irep-nburn,:)  = betas;
        sigma2draws(irep-nburn,:) = sigma2;
        etadraws(irep-nburn,:)    = eta;
        Vdraws(irep-nburn,:,:)    = V;        
        sdraws(irep-nburn,:)      = s';
        Pdraws(irep-nburn,:,:)    = P;
        
        %==========Forecasting
        s_fore = [s ; zeros(1,1)];
        betas_fore = [betas(1,(s(t)-1)*k+1:s(t)*k) ; zeros(1,k)];
        sigma2_fore = [sigma2(1,s(t)); zeros(1,1)];
        nii = sum(s_fore==s_fore(t));
        p_fore = betarnd(a+nii,b+1);    % Draw Pkk
        % Draw S(T+h)                         
        if p_fore > rand
            s_fore(t+1) = s_fore(t+1-1);       
        else
            s_fore(t+1) = s_fore(t+1-1) + 1;
        end
        
        % Update the parameters
        if s_fore(t+1) == s_fore(t)
            betas_fore(2,:) = betas_fore(1,:);
            sigma2_fore(2,:) = sigma2_fore(1,:);
        else
            betas_fore(2,:) = betas_fore(1,:) + randn(1,k)*chol(V);
            % We do not want to forecast with explosive
            % parameter draws. Here I check only the AR(1)
            % coefficient and if it is explosive I set it back
            % to the previous regime (instead of drawing again
            % until stationary which can be inefficient computationally)
            if sum(abs(betas_fore(2,1*constant + 1:end))) >= 1.000                      
                betas_fore(2,1*constant + 1) = betas_fore(1,1*constant + 1);
            end
            log_sigma_tplusone = log(sigma2_fore(1,:)) + sqrt(eta)*randn;
            sigma2_fore(2,:) = exp(log_sigma_tplusone);
        end

        y_hat = X_fore*betas_fore(2,:)' + sqrt(sigma2_fore(2))*randn;
        y_foredraws(irep-nburn,:) = y_hat;  
    end  % End saving after burn-in draws
end    %GIBBS iterations end here
fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8)
toc;