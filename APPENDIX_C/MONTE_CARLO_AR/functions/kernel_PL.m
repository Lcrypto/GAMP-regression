function [f] = kernel_PL(Z,X,bandw)
%Kernel density estimator
nobs = size(X,1);
n = size(Z,1); %#ok<NASGU>
Q1 = quantile(X,0.25);
Q3 = quantile(X,0.75);

if nargin==2
    h=0.79*(Q3-Q1)*nobs^(-0.2); 
else
    h = bandw;
end

 v1 = normpdf((Z - X)/h);
 f = mean(v1)/h;