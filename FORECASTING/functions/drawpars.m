function [Bdraw,vv,nu2,vvv] = drawpars(y,X,nu1,nu2,V,eta,t,k,M,s,dur) 	
%--------------------------------------------------------
% Draw the mean equation time-varying parameters Bt
% INPUTS:
% y  -  dependent variable 
% X  -  matrix of predictor variables
% nu -  measurement equation variance
% V  -  state equation covariance matrix
% t  -  time periods (sample)
% k  -  number of predictor variables in X
% s  -  regime indicators from Chib's algorithm
%--------------------------------------------------------

if nargin~=11
    error('Wrong number of inputs in drawpars')
end

t_m = zeros(t,k);  %storage matrix for state mean
t_v = zeros(k,t*k);  %storage matrix for state variance
c_m = zeros(k,1);   %initial state mean is 0
t_m(1,:) = c_m';
v_m = 10*eye(k);   %initial state variance is 100
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

