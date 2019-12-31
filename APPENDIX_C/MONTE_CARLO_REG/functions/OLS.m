function [beta,sigma] = OLS(X,Y)

beta = (X'*X)\(X'*Y);
sigma = (Y-X*beta)'*(Y-X*beta)/(size(Y,1)-size(X,2));