%% univariate stochastic volatility
function [h S] = SV(Ystar,h,S,phi,sig,mu,Dh)

T = length(S);

%% normal mixture
pi = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575]';
mi = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819]' - 1.2704;  %% means already adjusted!! %%
sigi = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261]';
sqrtsigi = sqrt(sigi);

%% sample S
temprand = rand(T,1);
for t = 1:T
    qi = pi.*normpdf(Ystar(t),h(t)+mi, sqrtsigi);
    qi = qi/sum(qi);
    S(t) = find(temprand(t)<cumsum(qi),1);
end
    
%% sample h
% y^* = h + d + \epison, \epison \sim N(0,\Omega),
% Hh = \alpha + \nu, \nu \ sim N(0,S),
% where d_t = Ez_t, \Omega = diag(\omega_1,\ldots,\omega_n), 
% \omega_t = var z_t, S = diag(Dh/(1-\phi^2), sig, \ldots, sig)
Hh =  speye(T,T) - spdiags(phi*ones(T-1,1),-1,T,T);
invSh = spdiags([(1-phi^2)/Dh; 1/sig*ones(T-1,1)],0,T,T);
dconst = mi(S); invOmega = spdiags(1./sigi(S),0,T,T);
h0 = Hh\[mu; ((1-phi)*mu)*ones(T-1,1)];
Kh = Hh'*invSh*Hh;
Ph = Kh + invOmega;
Ch = chol(Ph);
hhat = Ph\( Kh*h0 + invOmega*(Ystar-dconst));
h = hhat + Ch\randn(T,1);

