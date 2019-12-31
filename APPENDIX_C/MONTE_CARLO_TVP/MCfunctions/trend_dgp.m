function [y,x,theta_t,sigma] = trend_dgp(T)

% **************************************************************************************************************************
% Written by Dimitris Korobilis on 22/04/2017
% University of Essex
% **************************************************************************************************************************

%% Check for INPUT arguments
if nargin == 0
    T   = 200;            % Time series observations
end
mu = 4*rand;
sigma = 1;

%% Generate sparse time-varying parameter model
x=ones(T,1);
% 2. Generate unrestricted regression coefficients theta_t
% theta_t = zeros(T+100,1);
% for t = 1:T+100
%     if t == 1
%         theta_t(t,:) = mu + (1/(T^(2/4)))*randn;
%     else
%         theta_t(t,:) = theta_t(t-1,:) + (1/(T^(2/4)))*randn;
%     end
% end

% theta_t = zeros(T+100,1);
% z = 1*rand(T+100,10);
% gamma = 2*rand(10,1)-1;
% for t = 1:T+100
%     theta_t(t,:) = mu +  z(t,:)*gamma + (1/(T^(3/4)))*randn;
% end

for t = 1:T+100
    theta_t(t,:) = mu + (sign(2*rand-1)*(mu)*poissrnd(.1)) + (1/(T^(3/4)))*randn;
end

theta_t = theta_t(101:end,:);

% 4. Generate dependent variable y
y = sum(x.*theta_t,2) + sqrt(sigma).*randn(T,1);

