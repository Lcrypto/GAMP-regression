%% initialize TVD1
function [gam K h S Sig2 sig phi mu c0] = init1(W,X,Y,Srep)

rand('state', sum(100*clock) ); randn('state', sum(200*clock) );

[T m] = size(W);
p = size(X,2);
Z = [ones(60,1) W(1:60,:) X(1:60,:)];
ols0 = (Z'*Z)\(Z'*Y(1:60,:)); %% initial values
q = p+m+1;  % dim of states
k = p+2; % no. of vales S can take on
Tq = T*q;

sig0 = sum((Y(1:60,:) - Z*ols0).^2)/60;
h = log(sig0)*ones(T,1);
Sig2 = diag(.002*ones(q,1));    
invSig2 = Sig2\eye(q);
Dh = 1;  D = 5*speye(q); invD = inv(D);
H = speye(Tq,Tq) - [ [ sparse(q, (T-1)*q); kron(speye(T-1,T-1), speye(q))] sparse(Tq, q)];
gam0 = H\[ols0; zeros((T-1)*q,1)];
gam = gam0;
tempS = rand(T,1); S = zeros(T,1);
cumsumpi = cumsum([0.0073 .10556 .00002 .04395 .34001 .24566 .2575]);
for t=1:T
    S(t) = find(tempS(t)<cumsumpi,1);
end
phi1 = 20; phi2 = 1.5; 
phi = .8; sig = 1; mu = mean(h);
cbeta0 = 10; calp0 = 2/k*cbeta0/(1-2/k);
g = @(x) (phi1-1)*log((1+x)/2) + (phi2-1)*log((1-x)/2);
K = zeros(T,p); like = zeros(T,k); F = zeros(T,k); Kid = zeros(T,1);
c0 = 1/k; P = (1-c0)/(k-1)*ones(k,k) + (c0-(1-c0)/(k-1))*eye(k);
for loop=1:500
    %% sample K
    for t=1:T
        gamt = gam((t-1)*q+1:t*q);
        for j=1:k
            xt = [1 W(t,:) X(t,:).*Srep(j,:)]';
            u = Y(t)-xt'*gamt;
            like(t,j) = -.5*(log(2*pi)+h(t)) - .5*exp(-h(t)) * u^2;
        end       
    end
    like = exp(like)+1e-10;   
    F(1,:) = ones(1,k)/k; 
    for t=2:T
        Ft = (F(t-1,:)*P).*like(t,:);
        F(t,:) = Ft/sum(Ft);
    end
    Kid(T) = find(cumsum(F(T,:))>rand,1);
    K(T,:) = Srep(Kid(T),:);
    for t=T-1:-1:1
        prSt = F(t,:)' .* P(:,Kid(t+1));
        prSt = prSt/sum(prSt);
        Kid(t) = find(cumsum(prSt)>rand,1);
        K(t,:) = Srep(Kid(t),:);
    end   
    
    %% sample c0  % parameter for the Markov prior of K_t
    tempI = bin2dec(num2str(K));
    It = (tempI(2:end) == tempI(1:end-1));
    c0 = betarnd(sum(It)+calp0, T-1-sum(It)+cbeta0);
    P = (1-c0)/(k-1)*ones(k,k) + (c0-(1-c0)/(k-1))*eye(k);
    
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
    
    %% sample h and S
    Ystar = log((Y-bigX*gam).^2 + .0001 );
    [h S] = SV(Ystar,h,S,phi,sig,mu,Dh);
    
    %% sample Sig2 (diagonal matrix)
    err2 = reshape(H*gam,q,T);
    diaginvSig2 = gamrnd((T-1)/2, 2./sum(err2(:,2:end).^2,2));
    invSig2 = diag( diaginvSig2 );
    Sig2 = diag(1./diaginvSig2);    
    
    %% sample sig
    newSh = (h(1)-mu)^2*(1-phi^2) + sum(((h(2:T)-mu) - phi*(h(1:T-1)-mu)).^2);
    sig = 1/gamrnd(T/2, 2/newSh);
    
    %% sample phi
    tempsum = sum((h(1:T-1)-mu).^2);
    phihat = (h(2:T) - mu)'*(h(1:T-1)-mu)/tempsum;
    Vphi = sig / tempsum;
    phic = phihat + sqrt(Vphi)*randn;
    if abs(phic)<.999
        alp = exp(g(phic)-g(phi));
        if alp>rand
            phi = phic;          
        end
    end 
    
    %% sample mu
    Vmu = sig/((T-1)*(1-phi)^2 + (1-phi^2));
    muhat = Vmu * ((1-phi^2)/sig*h(1) + (1-phi)/sig*sum(h(2:end) - phi*h(1:T-1)));
    mu = muhat + sqrt(Vmu) * randn; 
end

