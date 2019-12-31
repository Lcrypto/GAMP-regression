function [y_fore] = TVP_TVD(y,x,x_fore,plag,nloop,burnin)

Y = y;
W = x(:,2:1+plag);
X = x(:,plag+2:end);
N = 10;

m = plag; % no of lags
T = size(Y,1);
p = size(X,2);
q = p+m+1;  % dim of states
k = p+2; % no. of vales S can take on
Tq = T*q; q2 = q^2;

%% initialize
countphi=0;
shpara = zeros(nloop - burnin,3);
sb0 = zeros(nloop-burnin,p);
sc0 = zeros(nloop-burnin,1);
sh = zeros(nloop-burnin,T);
sgam = zeros(nloop-burnin,Tq);
sK = zeros(nloop-burnin,T*p);
sgamtK = zeros(Tq,1); sgamtK2 = zeros(Tq,1);
sSig2 = zeros(nloop-burnin,q);
smut = zeros(T,q);
sOmegat = zeros(T,q,q);

temppri = zeros(k,1);
tempEYtp1 = zeros(N,1); tempprelike = zeros(N,1);

%% prior
balp0 = 1.5; bbeta0 = 1.5;
cbeta0 = 2; calp0 = 1;
nu02 = 1; S02 = [.002*ones(m+1,1); .002*ones(p,1)];  % for independent inverse-Gamma prior
nuh0 = 1; Sh0 = 0.01*(nuh0/2-1)*2;
phi1 = 20; phi2 = 1.5;
Dh = 5;  D = 5*speye(q); invD = inv(D);
H = speye(Tq,Tq) - [ [ sparse(q, (T-1)*q); kron(speye(T-1,T-1), speye(q))] sparse(Tq, q)];
% Z = [ones(60,1) W(1:60,:) X(1:60,:)];
% ols = (Z'*Z)\(Z'*Y(1:60)); %% prior
ols = zeros(q,1);
% ols = (Z'*Z)\(Z'*Y(1:60)); ols(1:m+1) = 0;
gam0 = H\[ols; zeros((T-1)*q,1)];

%% construct a few things
g = @(x) (phi1-1)*log((1+x)/2) + (phi2-1)*log((1-x)/2);
Srep = zeros(k,p);  % all the possible values of S
Srep(2:end-1,:) = eye(p);
Srep(end,:) = ones(1,p);

newnuh = nuh0 + T;
newnu2 = nu02 + T-1;

[gam,K,h,S,Sig2,sig,phi,mu,c0] = init1(W,X,Y,Srep);
invSig2 = Sig2\eye(q);

y_fore = zeros(nloop,1); % forecast
tic;
fprintf('Now you are running TVP_TVD')
fprintf('\n')
fprintf('Iteration 0000')
savedraw = 0;
for loop = 1:nloop+burnin
    if mod(loop,500) == 0
        fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%s%4d',8,8,8,8,8,8,8,8,8,8,8,8,8,8,'Iteration ',loop)
    end
    %% sample b0
    b0 = betarnd(K(1,:)+balp0, 1-K(1,:)+bbeta0)';
    
    %% sample c0  % parameter for the Markov prior of K_t
    tempI = bin2dec(num2str(K));
    It = (tempI(2:end) == tempI(1:end-1));
    c0 = betarnd(sum(It)+calp0, T-1-sum(It)+cbeta0);

    %% sample K
    K = GCK_reg1_markov(W,X,Y,K,Srep,h,Sig2,ols,D,b0,c0);
        
    %% sample gam
    bigX = SURform([ones(T,1) W X.*K]);
    invS = kron(speye(T,T), invSig2); invS(1:q,1:q) = invD;
    Q = H'*invS*H;
    XinvSig1 = bigX'*sparse(1:T,1:T,exp(-h));
    XinvSig1X = XinvSig1*bigX;
    invP = Q + XinvSig1X;
    Cp = chol(invP);
    gamhat = invP\(Q*gam0 + XinvSig1*Y);
    gam = gamhat + Cp\randn(Tq,1); 
    
    %% sample Sig2 (diagonal matrix)
    err2 = reshape(H*gam,q,T);
    newS2 = S02 + sum(err2(:,2:end).^2,2);
    diaginvSig2 = gamrnd(newnu2/2, 2./newS2);
    invSig2 = diag( diaginvSig2 );
    Sig2 = diag(1./diaginvSig2);
    
    % sample h
    Ystar = log((Y-bigX*gam).^2 + .0001 );
    h = SV(Ystar,h,S,phi,sig,mu,Dh);
    
    %% sample sig
    newSh = Sh0 + (h(1)-mu)^2*(1-phi^2) + sum(((h(2:T)-mu) - phi*(h(1:T-1)-mu)).^2);
    sig = 1/gamrnd(newnuh/2, 2/newSh);
    
    %% sample phi
    tempsum = sum((h(1:T-1)-mu).^2);
    phihat = (h(2:T) - mu)'*(h(1:T-1)-mu)/tempsum;
    Vphi = sig / tempsum;
    phic = phihat + sqrt(Vphi)*randn;
    if abs(phic)<.999
        alp = exp(g(phic)-g(phi));
        if alp>rand
            phi = phic;
            countphi = countphi+1;
        end
    end 
    
    %% sample mu
    Vmu = sig/((T-1)*(1-phi)^2 + (1-phi^2));
    muhat = Vmu * ((1-phi^2)/sig*h(1) + (1-phi)/sig*sum(h(2:end) - phi*h(1:T-1)));
    mu = muhat + sqrt(Vmu) * randn;
           
    if loop>burnin
        i=loop-burnin;
%         if mod(i,6)==0
%             savedraw = savedraw+1;
        sgam(i,:) = gam';
        sgamtK = sgamtK + gam.*reshape([ones(T,m+1) K]',T*q,1);
        sgamtK2 = sgamtK2 + (gam.*reshape([ones(T,m+1) K]',T*q,1)).^2;
        
        sK(i,:) =reshape(K',T*p,1);
        sh(i,:) = h';
        shpara(i,1) = phi;
        shpara(i,2) = sig;        
        shpara(i,3) = mu;
        sSig2(i,:) = diag(Sig2)';
        sb0(i,:) = b0';
        sc0(i,:) = c0;
        
        %% prediction - get N draws at a time
        htp1 = mu + phi*(h(end)-mu) + sqrt(sig)*randn(N,1);
        gamtp1 = repmat(gam(end-q+1:end),1,N) + repmat(sqrt(diag(Sig2)),1,N).*randn(q,N);
        sigtp1 = exp(htp1);
        P = (1-c0)/(k-1)*ones(k,k) + (c0-(1-c0)/(k-1))*eye(k);
        j = lookupKt(K(end,:));
        cumsumP = cumsum(P(j,:));
        for j=1:N
            Ktp1 = Srep(find(cumsumP > rand,1),:)';
            tempEYtp1(j) = x_fore*[gamtp1(1:m+1,j); Ktp1.*gamtp1(m+2:end,j)];
        end
        EYtp1 = mean(tempEYtp1);
        
        y_fore(i,:) = EYtp1; 
%         end
    end    
end
fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8)
% fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8)
toc; % Stop timer and print total time

