function P=drawPun(y,X,M,s,a,b)

% PURPOSE: computes one step of the Gibbs sampling for the transition prob
% matrix WITHOUT IMPOSING THAT THE LAST REGIME IS AN ABSORBING STATE!!!


[N,k]=size(X); k=k-1;
P=zeros(M+1,M+1);
%Compute the number of one-step transitions from state i to state i in
%the sequence s
clear nii;
for i=1:M+1
    nii(i)=sum(s==i);
    % Fill in the P matrix at the jth Gibbs sampling interation
    P(i,i) = betarnd(a+nii(i),b+1);
    if i<=M
        P(i,i+1)=1-P(i,i);  %fulfill the row restrictions
    end
end
