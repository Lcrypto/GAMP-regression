function [kdraw,lpy2n] = gck(yg,hh,capg,capf,sigv,kold,t,ex0,vx0,nvalk,kprior,kvals,p,kstate)
%GCK Function to implement Gerlach, Carter and Kohn (2000), 'Efficient 
%    Bayesian Inference for Dynamic Mixture Models', JASA. 
%--------------------------------------------------------------------------
% The dynamic mixture model is of the form
%    
%        Y[t] = h[t] x X[t] + G[t] x u[t]    
%        X[t] = F[t] x X[t-1] + GAMMA[t] x K[t] x v[t]
%
% where u[t] and v[t] are the errors (normal or conditionally normal).
% Conmditional on the parameters g[t], h[t], G[t], f[t], F[t] amd GAMMA[t],
% this algorithm provides estimates of the indicators K[t] as in Gerlach,
% Carter and Kohn (GCK), Recursion 1, pages 821-822.
%--------------------------------------------------------------------------
% Inputs: yg = p by t data matrix
%         hh = tp by m matrix (note that this is transpose of GKCs definition)
%         capg = t*p by p matrix containing gt which is akin to the st dev of measurement equation
%         capf = tm by m matrix from state equation (set to identies for random walk evolution of states)
%         sigv is the standard deviation (or sigv*sigv') is the state equation error variance when 
%                      regime change occurs (i.e. Kt=1). sigv will be an m by m matrix 
%         kold is the previous draw of k, which is t by 1 (i.e. this code only allows for one k)
%         t = nunber of observations
%         kstate = dimensionality of state vector
%         ex0, vx0 = mean and variance for initial values for state equation (m by 1 and m by m)
%         nvalk = number of values k can take on -- usually 2 for 0/1
%         kprior = prior probabilities for each value for k = nvalk by 1 (this is okay for 
%                      Bernoulli case, but in general may make this nvalk by t)       
%         kvals are the values k can take on -- usually 0/1           


% GCK's Step 1 on page 821
lpy2n=0;
mu = zeros(t*kstate,1);
omega = zeros(t*kstate,kstate);
for i = t-1:-1:1
    gatplus1 = sigv*kold(i+1,1);
    ftplus1 = capf(kstate*i+1:kstate*(i+1),:);
    cgtplus1 = capg(i*p+1:(i+1)*p,:);
    htplus1 = hh(i*p+1:(i+1)*p,:)';
    rtplus1 = (htplus1'*gatplus1)*(htplus1'*gatplus1)' + cgtplus1*cgtplus1';
    hrtinv = htplus1/(rtplus1);
    btplus1 = gatplus1*gatplus1'*hrtinv;
    atplus1 = (eye(kstate) - btplus1*htplus1')*ftplus1;
    if kold(i+1,1) == 0;
        ctplus1 = zeros(kstate,kstate);
    else
        cct = gatplus1*(eye(kstate) - gatplus1'*hrtinv*htplus1'*gatplus1)*gatplus1';
        ctplus1 = chol(cct)';
    end
    otplus1 = omega(kstate*i+1:kstate*(i+1),:);
    dtplus1 = ctplus1'*otplus1*ctplus1 + eye(kstate);
    ocdc = otplus1*(ctplus1/dtplus1)*ctplus1';
    omega(kstate*(i-1)+1:kstate*i,:) = atplus1'*(otplus1 - ocdc*otplus1)*atplus1 + ftplus1'*hrtinv*htplus1'*ftplus1;
    mutplus1 = mu(kstate*i+1:kstate*(i+1),:);
    mu(kstate*(i-1)+1:kstate*i,:) = atplus1'*(eye(kstate) - ocdc)*(mutplus1 - otplus1*btplus1*yg(:,i+1)) + ftplus1'*hrtinv*yg(:,i+1);  
end

% GCKs Step 2 on pages 821-822
kdraw = kold;
ht = hh(1:p,:)';
ft = capf(1:kstate,:);
gat = zeros(kstate,kstate);
% Note: this specification implies no shift in first period -- sensible
rt = ht'*ft*vx0*ft'*ht + ht'*gat*gat'*ht+ capg(1:p,:)*capg(1:p,:)';
rtinv = inv(rt);
jt = (ft*vx0*ft'*ht + gat*gat'*ht)*rtinv;
mtm1 = (eye(kstate) - jt*ht')*(ft*ex0) + jt*(yg(:,1));
vtm1 = ft*vx0*ft'+ gat*gat' - jt*rt*jt';
lprob = zeros(nvalk,1);
for i = 2:t
    ht = hh((i-1)*p+1:i*p,:)';
    ft = capf(kstate*(i-1)+1:kstate*i,:);
    for j = 1:nvalk
        gat = kvals(j,1)*sigv;
        ggg = gat*gat';
        rt = ht'*ft*vtm1*ft'*ht + ht'*ggg*ht + capg((i-1)*p+1:i*p,:)*capg((i-1)*p+1:i*p,:)';
        rtinv = inv(rt);
        jt = (ft*vtm1*ft'*ht + ggg*ht)*rtinv;
        mt = (eye(kstate) - jt*ht')*(ft*mtm1) + jt*(yg(:,i));
        vt = ft*vtm1*ft' + ggg - jt*rt*jt';
        yghtg = yg(:,i) - ht'*(ft*mtm1);
        lpyt = -.5*log(det(rt)) - .5*(yghtg)'*rtinv*(yghtg);
        if det(vt)<=0;
            tt = zeros(kstate,kstate);
        else
            tt = chol(vt)';
        end
        ot = omega(kstate*(i-1)+1:kstate*i,:);
        mut = mu(kstate*(i-1)+1:kstate*i,:);
        tempv = eye(kstate) + tt'*ot*tt;
        lpyt1n = -.5*log(det(tempv)) -.5*(mt'*ot*mt - 2*mut'*mt - (mut - ot*mt)'*(tt/(tempv))*tt'*(mut - ot*mt));
        lprob(j,1) = log(kprior(j,1)) + lpyt1n + lpyt;
        if i == 2;
            lpy2n = lpyt1n + lpyt;
        end
    end
    pprob = exp(lprob)./sum(exp(lprob));
    tempv = rand(1,1);
    tempu = 0;
    for j = 1:nvalk
        tempu = tempu + pprob(j,1);
        if tempu > tempv;
            kdraw(i,1) = kvals(j,1);
            break
        end
    end

    gat = kdraw(i,1)*sigv;
    ggg = gat*gat';
    rt = ht'*ft*vtm1*ft'*ht + ht'*ggg*ht + capg((i-1)*p+1:i*p,:)*capg((i-1)*p+1:i*p,:)';
    rtinv = inv(rt);
    jt = (ft*vtm1*ft'*ht + ggg*ht)*rtinv;
    mtm1 = (eye(kstate) - jt*ht')*(ft*mtm1) + jt*(yg(:,i));
    vtm1 = ft*vtm1*ft' + ggg - jt*rt*jt';           
end