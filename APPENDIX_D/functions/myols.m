
function [bhat, seb] = myols(Y, X)
% Simple OLS regression stats
T = size(X,1);
K = size(X,2);

bhat = inv(X'*X)*(X'*Y);          % Coefficients  (also (X'*X)^-1*(X'*Y) or X\Y)
yhat = X*bhat;                    % XB or beta(1)*x1 + beta(2)*x2
u = Y - yhat;                     % Residuals
RSS = u'*u;                       % Sum of Squared Residuals 
s2 = RSS/(T-K);                   % Estimate of variance of residuals
s  = s2^0.5;                      % Standard error of estimate
covb = s2*inv(X'*X);              % Covariance of beta
seb = sqrt(diag(covb));           % se of beta
tb = bhat ./ seb;                 % t-stats of beta
pb = 1 - cdf('t', abs(tb), T-K);  % p-values
pb = 2*pb;                        % 2-tailed p-values

ym = Y - mean(Y);
TSS = ym'*ym;
R2 = 1.0 - RSS/TSS;                    % R-squared
R2adj = 1.0 - (RSS/(T-K))/(TSS/(T-1)); % Adj. R-squared

ESS = TSS - RSS;
F = (ESS/(K-1))/(RSS/(T-K));           % F-stat
pF = 1 - cdf('f', F, K-1, T-K);        % p-value of F