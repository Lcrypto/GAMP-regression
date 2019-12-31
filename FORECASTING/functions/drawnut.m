function [Sigtdraw,vv] = drawnut(y,X,vart,eta,t,k,s) 	
%--------------------------------------------------------
% Draw the mean equation time-varying parameters Bt
% INPUTS:
% y  -  dependent variable 
% X  -  matrix of predictor variables
% nu -  measurement equation variance
% eta  -  state equation covariance matrix
% t  -  time periods (sample)
% k  -  number of predictor variables in X
% s  -  regime indicators from Chib's algorithm
%--------------------------------------------------------

if nargin~=7
    error('Wrong number of inputs in drawnut')
end

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
        pst_vr = v_m + eta; %Predicted state covariance       
    else
        %no switch - do not update predicted covariance with sigma
        pst_vr = v_m; %Predicted state covariance
    end
    % Kalman filter iterations
    f_cast = y(i,:) - h*pst_cf;  %cond. forecast error
    ss = h*pst_vr*h' + vart(i,1);
    K_GN = pst_vr*h'*inv(ss); %Kalman gain
    c_m = pst_cf + K_GN*f_cast;  %updated state mean
    v_m = (eye(k) - K_GN*h)*pst_vr;   %updated state variance
    t_m(i,:) = c_m';     %store state mean
    t_v(:,(i-1)*k+1:i*k) = v_m;  %store state variance  
end

tt_m = t_m;
tt_v = t_v;
% Draw the last (at time T) theta, using mean and variance from forward iterations above

Sigtdraw(s(t),:) = tt_m(t,:) + sqrt(tt_v(:,t))*randn; %#ok<*AGROW>
   
if s(t,:)<t;       
    for i = s(t)+1:t
        Sigtdraw(i,:)=Sigtdraw(i-1,:) + sqrt(eta)*randn;
        vv((i-2)*k+1:(i-1)*k,:)=(Sigtdraw(i,:)-Sigtdraw(i-1,:))'*(Sigtdraw(i,:)-Sigtdraw(i-1,:));		             
    end
end

%start backward iterations
for i = t:-1:2       
    if s(i) ~= s(i-1)
        f_cast0m = Sigtdraw(s(i),:) - tt_m(i-1,:);
        ss0m = tt_v(:,i-1) + eta;
        K_GN0m = tt_v(:,i-1)*inv(ss0m);
        tt_m(i-1,:) = tt_m(i-1,:) + f_cast0m*K_GN0m;        
        tt_v(:,i-1) = (eye(k)-K_GN0m)*tt_v(:,i-1);
        Sigtdraw(s(i-1),:) = tt_m(i-1,:) + sqrt(tt_v(:,i-1))*randn;
        vv((s(i-1)-1)*k+1:s(i-1)*k,:) = (Sigtdraw(s(i),:)-Sigtdraw(s(i-1),:))'*(Sigtdraw(s(i),:)-Sigtdraw(s(i-1),:));      
    end
end
