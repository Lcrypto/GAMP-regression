function [Z] = orthogonalize(X,constant)


r = cholcov(cov(X(:,1+constant:end))); 
W = eye(size(X,2)); W(1+constant:end,1+constant:end)=r; 
Z = X/W; Z(:,1+constant:end) = Z(:,1+constant:end) - mean(Z(:,1+constant:end));