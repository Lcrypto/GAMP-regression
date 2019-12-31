function [pst_n,s_n,llikf]=drawstates(y,X,M,betas_o,sigma2_o,P_o)
 
% PURPOSE: computes one step of the Gibbs sampling for the states

[N,k]=size(X); k=k-1;

%Some initializations
pst_n=zeros(M+1,N); %prob of the states draw
s_n=zeros(N,1);       % states draw
pst_ytl=zeros(M+1,N);
pst_yt=zeros(M+1,N);

% 1. Calculate and store for t=1,...,n  p(st|yt,M,theta,P)
ilikf=zeros(N,M+1);   %store the individual likelihood function under all regimes

%Compute p(s1|Y1,theta,P) first

pst_ytl(:,1)=eye(M+1,1); %degenerate prob distr for p(s1|Y0,theta,P)
fyt_ytl=normgen(y(1),X(1,:),betas_o,sigma2_o,M); 
ilikf(1,:)=fyt_ytl';
num=(pst_ytl(:,1).*fyt_ytl);
den=sumc(num);
pst_yt(:,1)=num./den;

%Now repeat the same step for t=2,...,N
for t=2:N
    pst_ytl(:,t)=P_o'*pst_yt(:,t-1);
    %compute f(yt|Yt-1,theta_k) for each state k=1,...,M+1
    fyt_ytl=normgen(y(t),X(t,:),betas_o,sigma2_o,M); 
    ilikf(t,:)=fyt_ytl';
	num=(pst_ytl(:,t).*fyt_ytl);
    den=sumc(num);
    pst_yt(:,t)=num./den;
end

% Here I compute the likelihood function for t=1,..,n at each iteration j

likf = diag(ilikf*pst_ytl);
llikf = sum(log(likf));

% 2. Sample sn from p(sn|y,M,theta,P)
pst_n(:,N) = flipud(eye(M+1,1)); %degenerate at sn=M+1

% 3. Sample for t=n-1,n-2,...,1 st from p(st|y,M,S^t+1,theta,P)

for t=N-1:-1:2
    s_n(t+1)=discrete(pst_n(:,t+1));     % procedure to generate a discrete random variable
    num=pst_yt(:,t) .* P_o(:,s_n(t+1));
    pst_n(:,t)=num ./ sum(num);
end
s_n(t)=discrete(pst_n(:,t));

pst_n(:,1)=eye(M+1,1); %degenerate at s1=1
s_n(1)=1;  %at time 1 the state is equal to 1
