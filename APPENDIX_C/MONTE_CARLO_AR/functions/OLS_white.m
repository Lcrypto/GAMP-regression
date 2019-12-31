
function [estimator,SE,SE_robust,sigma2,resid] = OLS_white(Y,X,add_constant)

% Purpose: 
% Ordinary least squares with White robust standard error
% -----------------------------------
% Model:
% Yi = Xi * Beta + ui , where ui ~ N(0,s^2)
% -----------------------------------
% Algorithm: 
% inv(X'*X)* (X'*Y)
% -----------------------------------
% Usage:
% Y = dependent variable (n * 1 vector)
% X = regressors (n * k matrix)
% add_constant = whether to add a constant to X (default = 0)
% -----------------------------------
% Returns:
% estimator = estimator corresponding to the k regressors
% SE = standard error of the estimator
% SE_robust = White robust standard error of the estimator
% sigma2 = estimated variance of disturbances
% resid = residuals series of the regression
% In the absence of returning arguments, estimation results will be displayed on screen
% 
% Written by Hang Qian, Iowa State University
% Contact me:  matlabist@gmail.com



% Written by Hang Qian, Iowa State University

if nargin < 2;    error('Incomplete data');end
if nargin < 3;    add_constant = 0;end

[nrow_x,ncol_x] = size(X);
[nrow_y,ncol_y] = size(Y);
if nrow_x < ncol_x;    X = X';    ncol_x = nrow_x;end
if nrow_y < ncol_y;    Y = Y';    ncol_y = nrow_y;end
if ncol_x < ncol_y;    Y_temp = Y;    Y = X;    X = Y_temp;end


[nobs, nvar] = size(X);
if add_constant == 1
    disp('A constant is added to X')
    X = [ones(nobs,1),X];
    nvar = nvar + 1;
end


XX = X'*X;
estimator = XX \ (X'*Y);
resid = Y - X * estimator;
RSS = resid' * resid;
sigma2 = RSS / (nobs-nvar-1);
inv_XX = inv(XX);
cov_mat = inv_XX * sigma2;
SE = sqrt(diag(cov_mat));

X_trans = X .* resid(:,ones(1,nvar));
cov_mat_robust = inv_XX * (X_trans' * X_trans) * inv_XX;
cov_mat_robust = nobs / (nobs-nvar-1) * cov_mat_robust; 
SE_robust = sqrt(diag(cov_mat_robust));

t_stat = estimator ./ SE_robust;
try
    p_val = (1 - tcdf(abs(t_stat),nobs-nvar-1)) * 2;
catch
    p_val = (1 - 0.5 * erfc(-0.7071 * abs(t_stat))) * 2;
end

Y_demean = Y - mean(Y);
R2 = 1-RSS/(Y_demean' * Y_demean);
log_like = -nobs/2 * (1 + log(2*pi) + log(RSS/nobs));
DW = sum(diff(resid).^2) / RSS;
eval([char([81 72 49 61]),'[87 114 105 116 116 101 110 32 98 121];'])
eval([char([81 72 50 61]),'[32 72 97 110 103 32 81 105 97 110];'])

if nargout == 0
    result = cell(nvar+1,6);
    result(1,:)={'',' Estimator','SE','SE_robust','t-stat','p-val'};
    for m=1:nvar
           result(m+1,1)={['C(',num2str(m-add_constant),')']};
           result(m+1,2:end)={estimator(m),SE(m),SE_robust(m),t_stat(m),p_val(m)};
    end
    disp(' ')
    disp(['Number of observations ',num2str(nobs)])
    disp(['Number of regressors ',num2str(nvar)])
    disp(result)
        
    disp(['R2:    ',num2str(R2)])
    disp(['Variance of s2:    ',num2str(sigma2)])
    disp(['Log likelihood: ',num2str(log_like)])
    disp(['D-W statistics:      ',num2str(DW)])
    fwrite(1, char([QH1,QH2,10,13]))
end


end

