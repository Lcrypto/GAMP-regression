function [beta,sigma2,bcov,t_stats] = ols_quantities(y,x,intercept)

[n,k] = size(x);

if intercept == 0
elseif intercept == 1
    x = [ones(n,1) x];
else
    error('Wrong specification of intercept. Please select 0/1 only.')    
end

beta = inv(x'*x)*(x'*y);       % OLS estimate of bean
e = y - x*beta;             % OLS error
sse = e'*e;                    % OLS sum of squared errors
sigma2 = sse/(n-k-1);              % OLS estimate of variance
bcov = sigma2*inv(x'*x);           % OLS estimate of covariance of beta
bcovt = sqrt(diag(bcov));
t_stats = beta(2:end)./bcovt(2:end);



