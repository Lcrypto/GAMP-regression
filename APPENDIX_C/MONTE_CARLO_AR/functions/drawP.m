function P=drawP(y,X,M,s)

% PURPOSE: computes one step of the Gibbs sampling for the transition prob
% matrix

global a b;

[N,k]=size(X); k=k-1;
P=zeros(M+1,M+1);
%Compute the number of one-step transitions from state i to state i in
%the sequence s
clear nii;
for i=1:M
    nii(i)=sum(s==i);
    % Fill in the P matrix at the jth Gibbs sampling interation
    P(i,i)=betarnd(a+nii(i),b+1);
    P(i,i+1)=1-P(i,i);  %fulfill the row restrictions
end
nii(M+1)=N - nii*ones(M,1);
P(M+1,M+1)=1;  %prob of staying in states M+1 once reached equal to 1
