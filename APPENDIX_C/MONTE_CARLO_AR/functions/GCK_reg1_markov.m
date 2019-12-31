%% GCK algorithm for the first formulation of the TVD regression example 

function K = GCK_reg1_markov(W,X,Y,K,Srep,h,Sig2,ols,D,b0,c0)
% b0 is a kx1 vector of paraemters for the initial distn of K_1
% c0 controls the smoothness of transition from K_t to K_{t+1}

[T p] = size(K);
m = size(W,2);
q = p+m+1;  % dim of states
k = p+2; % no. of vales S can take on
smut = zeros(T,q);
sOmegat = zeros(T,q,q);
Omegatp1 = zeros(q,q); mutp1 = zeros(q,1);
pK = zeros(k,1);
for t=T-1:-1:1
    xtp1 = [1 W(t+1,:) X(t+1,:).*K(t+1,:) ]';        
    rtp1 = xtp1'*Sig2*xtp1 + exp(h(t+1));
    Btp1 = Sig2 * xtp1 / rtp1;
    Atp1 = speye(q) - Btp1*xtp1';
    Ctp1 = chol(Sig2 - Sig2*(xtp1*xtp1')*Sig2/rtp1)';
    Dtp1 = Ctp1' * Omegatp1 * Ctp1 + speye(q);
    invDtp1 = inv(Dtp1);
    sOmegat(t,:,:) = Atp1'*(Omegatp1-Omegatp1*Ctp1*invDtp1*Ctp1'*Omegatp1)*Atp1 +...
        xtp1*xtp1'/rtp1;
    smut(t,:) = Atp1'*(speye(q) - Omegatp1*Ctp1*invDtp1*Ctp1')*(mutp1-Omegatp1*Btp1*Y(t+1))+...
        xtp1*Y(t+1)/rtp1;
       
    Omegatp1 = squeeze(sOmegat(t,:,:));
    mutp1 =  smut(t,:)';
end
mtm1 = ols; Vtm1 = D;
for t=1:T
    Omegat = squeeze(sOmegat(t,:,:));
    mut = smut(t,:)';
    for j=1:k 
        xt = [1 W(t,:) X(t,:).*Srep(j,:)]';
        Rt = xt'*Vtm1*xt + xt'*Sig2*xt + exp(h(t));
        Jt = (Vtm1 +Sig2)*xt/Rt;
        mt = (speye(q) - Jt*xt')*mtm1 + Jt*Y(t);
        Vt = Vtm1 + Sig2 - Jt*Jt'*Rt;
        Tt = chol(Vt)';
        TOT = Tt'*Omegat*Tt + speye(q);
        TmutmOmt = Tt'*(mut-Omegat*mt);
            
        pK(j) = -log(Rt)/2 - (Y(t)-xt'*mtm1)^2/(2*Rt) + ...
            -log(det(TOT))/2 - 1/2*(mt'*Omegat*mt - 2*mut'*mt - TmutmOmt'* (TOT\TmutmOmt));
    end     
        
    %% Markov prior 
    if t==1
        ii = lookupKt(K(2,:));
        tempi = sum(repmat(log(b0'),k,1).*Srep + repmat(log(1-b0)',k,1).*(1-Srep),2) ;
        tempj = log((1-c0)/(k-1))*ones(k,1);  tempj(ii) = log(c0);
        logKpri = tempi + tempj;
    elseif t==T
        ii = lookupKt(K(T-1,:));
        logKpri = log((1-c0)/(k-1))*ones(k,1);
        logKpri(ii) = log(c0);
    else
        ii = lookupKt(K(t-1,:));
        j = lookupKt(K(t+1,:));
        tempi = log((1-c0)/(k-1))*ones(k,1); tempi(ii) = log(c0);
        tempj = log((1-c0)/(k-1))*ones(k,1); tempj(j) = log(c0);
        logKpri = tempi + tempj;
    end    
    pK1 = pK + logKpri;
    pK2 = exp(pK1-max(pK1));
    pK3 = pK2 /sum(pK2);
    ind = find(cumsum(pK3)>rand, 1 );
    K(t,:) = Srep(ind,:);

    %% update the parameters
    xt = [1 W(t,:) X(t,:).*K(t,:)]';
    Rt = xt'*Vtm1*xt + xt'*Sig2*xt + exp(h(t));
    Jt = (Vtm1 +Sig2)*xt/Rt;
    mtm1 = (speye(q) - Jt*xt')*mtm1 + Jt*Y(t);
    Vtm1 = Vtm1 + Sig2 - Jt*Jt'*Rt;    
 end   