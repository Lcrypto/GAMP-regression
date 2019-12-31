function [y,beta_t,V_t] = tvp_reg_dgp(T)

%% Generate sparse time-varying parameter model

% % 2. Assume theta_t has three levels, theta_1, theta_2, theta_3
% theta_1 = 3.5; theta_2 = -2.5; theta_3 = 2;
% theta_t = zeros(T,1);
% for t = 1:T
%     if t < round(T/3)
%         theta_t(t,:) = theta_1;
%     elseif t >= round(T/3) && t < 2*round(T/3)
%         theta_t(t,:) = theta_2;
%     else
%         theta_t(t,:) = theta_3;
%     end
% end
% beta_t = theta_t;

% 3. Generate unrestricted regression coefficients sigma_t
sigma_t = zeros(T+100,1);
sigma0 = log(1);  % Regression variance
for t = 1:T+100
    if t == 1
        sigma_t(t,:) = sigma0 + (1/(T^(2/4)))*randn;
    else
        sigma_t(t,:) = 0.99*sigma_t(t-1,:)  + (1/(T^(2/4)))*randn;
    end
end
V_t = exp(sigma_t(101:end,:));

% 4. Generate dependent variable y
y = beta_t + sqrt(V_t).*randn(T,1);

