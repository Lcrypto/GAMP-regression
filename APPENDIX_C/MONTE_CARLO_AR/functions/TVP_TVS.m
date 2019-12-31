function [beta,sigma] = TVP_TVS(target, data, mumean1b, mulambdastar, numbofits, burnin)

warning off all
target = target';
every = 1;
[T, p]=size(data);

%mumean1b = 2;
%mulambdastar = 0.1;

rho = 0.97 * ones(1, p);
rhobeta = 0.97 * ones(1, p);
lambda = mulambdastar * ones(1, p);
mu = mumean1b * ones(1, p);
Delta = rho ./ (1 - rho) .* lambda ./ mu;

start_samples = 500;
start_adap = 1000;

newbeta = zeros(T, p);
Psi = ones(T, p);
kappa = zeros(T, p);
for j = 1:p
    kappa(:, j) = poissrnd(Delta(j) * Psi(:, j));
end
kappasigmasq = ones(T, 1);
lambdasigma = 3;
musigma = 0.03;
rhosigma = 0.95;

sigmasq = musigma * ones(T, 1);

xtx = zeros(p, p, T);
for i = 1:T
    xtx(:, :, i) = data(i, :)' * data(i, :);
end

lambdastar = 1;
mustar = 1;

Psisd = 0.01 * ones(T, p);
loglambdasd = log(0.1) * ones(1, p);
logmeansd = log(0.1) * ones(1, p);
logscale = 0.5 * log(2.4^2 / 4) * ones(1, p);
logrhobetasd = log(0.1) * ones(1, p);
logrhosd = log(0.1) * ones(1, p);
logrhosigmasd = log(0.01);
loglambdasigmasd = log(0.001);
loggammasigmasqsd = log(0.001);
logscalesigmasq = 0.5 * log(2.4^2 / 3);
mugammasd = 0.003;
vgammasd = 0.003;
logsigmasqsd = log(0.001) * ones(1, T);
logkappaq = 4 / 3 * ones(T - 1, p);
logkappasigmasqq = 4 / 3 * ones(T, 1);

kappaaccept = 0;
kappacount = 0;
kappalambdasigmaccept = 0;
kappasigmasqcount = 0;
sigmasqparamaccept = zeros(1, p);
sigmasqparamcount = zeros(1, p);
Psiparamaccept = zeros(1, p);
Psiparamcount = zeros(1, p);
mugammaaccept = 0;
mugammacount = 0;
vgammaaccept = 0;
vgammacount = 0;
sigmasq1accept = zeros(1, T);
sigmasq1count = zeros(1, T);
numberofiterations = burnin + every * numbofits;

holdPsi = zeros(T, p, numbofits);
holdbeta = zeros(T, p, numbofits);
holdsigmasq = zeros(T, numbofits);
holdlambda = zeros(p, numbofits);
holdmu = zeros(p, numbofits);
holdrhobeta = zeros(p, numbofits);
holdrho = zeros(p, numbofits);
holdlambdasigma = zeros(1, numbofits);
holdmusigma = zeros(1, numbofits);
holdrhosigma = zeros(1, numbofits);
holdlambdastar = zeros(1, numbofits);
holdmustar = zeros(1, numbofits);

sum1 = zeros(4, p);
sum2 = zeros(10, p);
sum1sigmasq = zeros(3, 1);
sum2sigmasq = zeros(6, 1);

limit1 = 0.9999;

[meankf, varkf, loglike] = kf_NGAR(data, target, Psi, sigmasq, rhobeta);

cholstar = chol(varkf(:, :, T))';
newbeta(T, :) = (meankf(:, T) + cholstar * randn(size(cholstar, 2), 1))';
for i = (T-1):-1:1
    Gkal = diag(rhobeta .* sqrt(Psi(i + 1, :) ./ Psi(i, :)));
    invQ = diag(1 ./ (1 - rhobeta.^2) .* 1./Psi(i + 1, :));
    varfb = inv(inv(varkf(:, :, i)) + Gkal' * invQ * Gkal);
    meanfb = varfb * (inv(varkf(:, :, i)) * meankf(:, i) + Gkal' * invQ * newbeta(i + 1, :)');
    cholstar = chol(varfb)';
    newbeta(i, :) = (meanfb + cholstar * randn(size(cholstar, 2), 1))';
end

beta = newbeta(:, 1:p);

checkstar = 1;

tic;
fprintf('Now you are running TVP_TVS')
fprintf('\n')
fprintf('Iteration 0000')
for it = 1:numberofiterations
    if mod(it,500) == 0
        fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%s%4d',8,8,8,8,8,8,8,8,8,8,8,8,8,8,'Iteration ',it)
    end

%     if ( mod(it, 20) == 0 )
%         disp(['it = ' num2str(it)]);
%         disp(['lambdastar = ' num2str(lambdastar)]);
%         disp(['mustar = ' num2str(mustar)]);       
%         disp(['lambda sigmasq = ' num2str(lambdasigma)]);
%         disp(['mean sigmasq = ' num2str(musigma)]);
%         disp(['rhosigma = ' num2str(rhosigma)]);
%         disp(['lambda = ' num2str(lambda)]);
%         disp(['mu = ' num2str(mu)]);
%         disp(['rho beta = ' num2str(rhobeta)]);
%         disp(['rho = ' num2str(rho)]);
%         disp(' ')
%         disp(['kappa accept = ' num2str(kappaaccept/kappacount)]);
%         disp(['kappasigmasq accept = ' num2str(kappalambdasigmaccept/kappasigmasqcount)]);
%         disp(['Psi param accept = (' num2str(min(Psiparamaccept./Psiparamcount)) ', ' num2str(max(Psiparamaccept./Psiparamcount)) ')']);
%         disp(['sigmasq param accept = ' num2str(sigmasqparamaccept/sigmasqparamcount)]);
%         disp(['mugamma accept = ' num2str(mugammaaccept/mugammacount)]);
%         disp(['vgamma accept = ' num2str(vgammaaccept/vgammacount)]);
%         disp(['sigmasq accept = (' num2str(min(sigmasq1accept./sigmasq1count)) ', ' num2str(max(sigmasq1accept./sigmasq1count)) ')']);
%         disp(' ')
%     end
    
    
    
    if ( checkstar == 1 )
        
        for i = 1:T
            for j = 1:p
                newPsi = Psi(i, j) * exp(Psisd(i, j) * randn);
                
                if (i == 1 )
                    loglike = (lambda(j) - 1) * log(Psi(1, j)) - lambda(j) * Psi(1, j) / mu(j);
                    pnmean = Psi(i, j) * Delta(j);
                    loglike = loglike - pnmean + kappa(i, j) * log(pnmean);
                    loglike = loglike - 0.5 * log(Psi(1, j)) - 0.5 * beta(1, j)^2 / Psi(1, j);
                    var1 = Psi(2, j) * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(Psi(2, j) / Psi(1, j)) * beta(1, j);
                    loglike = loglike - 0.5 * (beta(2, j) - mean1)^2 / var1;
                    
                    newloglike = (lambda(j) - 1) * log(newPsi) - lambda(j) * newPsi / mu(j);
                    pnmean = newPsi * Delta(j);
                    newloglike = newloglike - pnmean + kappa(i, j) * log(pnmean);
                    newloglike = newloglike - 0.5 * log(newPsi) - 0.5 * beta(1, j)^2 / newPsi;
                    var1 = Psi(2, j) * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(Psi(2, j) / newPsi) * beta(1, j);
                    newloglike = newloglike - 0.5 * (beta(2, j) - mean1)^2 / var1;
                elseif ( i < T )
                    lam1 = lambda(j) + kappa(i - 1, j);
                    gam1 = lambda(j) / mu(j) + Delta(j);
                    loglike = (lam1 - 1) * log(Psi(i, j)) - gam1 * Psi(i, j);
                    pnmean = Psi(i, j) * Delta(j);
                    loglike = loglike - pnmean + kappa(i, j) * log(pnmean);
                    var1 = Psi(i, j) * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(Psi(i, j) / Psi(i - 1, j)) * beta(i - 1, j);
                    loglike = loglike - 0.5 * log(var1) - 0.5 * (beta(i, j) - mean1)^2 / var1;
                    var1 = Psi(i + 1, j) * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(Psi(i + 1, j) / Psi(i, j)) * beta(i, j);
                    loglike = loglike - 0.5 * (beta(i + 1, j) - mean1)^2 / var1;
                    
                    lam1 = lambda(j) + kappa(i - 1, j);
                    gam1 = lambda(j) / mu(j) + Delta(j);
                    newloglike = (lam1 - 1) * log(newPsi) - gam1 * newPsi;
                    pnmean = newPsi * Delta(j);
                    newloglike = newloglike - pnmean + kappa(i, j) * log(pnmean);
                    var1 = newPsi * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(newPsi / Psi(i - 1, j)) * beta(i - 1, j);
                    newloglike = newloglike - 0.5 * log(var1) - 0.5 * (beta(i, j) - mean1)^2 / var1;
                    var1 = Psi(i + 1, j) * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(Psi(i + 1, j) / newPsi) * beta(i, j);
                    newloglike = newloglike - 0.5 * (beta(i + 1, j) - mean1)^2 / var1;
                else
                    lam1 = lambda(j) + kappa(i - 1, j);
                    gam1 = lambda(j) / mu(j) + Delta(j);
                    loglike = (lam1 - 1) * log(Psi(i, j)) - gam1 * Psi(i, j);
                    var1 = Psi(i, j) * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(Psi(i, j) / Psi(i - 1, j)) * beta(i - 1, j);
                    loglike = loglike - 0.5 * log(var1) - 0.5 * (beta(i, j) - mean1)^2 / var1;
                    
                    lam1 = lambda(j) + kappa(i - 1, j);
                    gam1 = lambda(j) / mu(j) + Delta(j);
                    newloglike = (lam1 - 1) * log(newPsi) - gam1 * newPsi;
                    var1 = newPsi * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(newPsi / Psi(i - 1, j)) * beta(i - 1, j);
                    newloglike = newloglike - 0.5 * log(var1) - 0.5 * (beta(i, j) - mean1)^2 / var1;
                end
                
                
                logaccept = newloglike - loglike + log(newPsi) - log(Psi(i, j));
                accept = 1;
                if ( (isnan(logaccept) == 1) || (isinf(logaccept) == 1) )
                    accept = 0;
                elseif ( logaccept < 0 )
                    accept = exp(logaccept);
                end
                
                Psisd(i, j) = Psisd(i, j) + (accept - 0.3) / it^0.6;
                
                if ( rand < accept )
                    Psi(i, j) = newPsi;
                end
            end
        end
        
        
        
        for i = 1:(T-1)
            for j = 1:p
                newkappa = kappa(i, j) + (2 * ceil(rand) - 1) * geornd(1 / (1 + exp(logkappaq(i, j))));

                if ( newkappa < 0 )
                    accept = 0;
                else
                    lam1 = lambda(j) + kappa(i, j);
                    gam1 = lambda(j) / mu(j) + Delta(j);
                    loglike = lam1 * log(gam1) - gammaln(lam1) + (lam1 - 1) * log(Psi(i + 1, j));
                    pnmean = Psi(i, j) * Delta(j);
                    loglike = loglike + kappa(i, j) * log(pnmean) - gammaln(kappa(i, j) + 1);
                    
                    lam1 = lambda(j) + newkappa;
                    gam1 = lambda(j) / mu(j) + Delta(j);
                    newloglike = lam1 * log(gam1) - gammaln(lam1) + (lam1 - 1) * log(Psi(i + 1, j));
                    pnmean = Psi(i, j) * Delta(j);
                    newloglike = newloglike + newkappa * log(pnmean) - gammaln(newkappa + 1);
                    
                    logaccept = newloglike - loglike;
                    
                    accept = 1;
                    if ( (isnan(logaccept) == 1) || (isinf(logaccept) == 1) )
                        accept = 0;
                    elseif ( logaccept < 0 )
                        accept = exp(logaccept);
                    end
                end
                
                kappaaccept = kappaaccept + accept;
                kappacount = kappacount + 1;
                
                if ( rand < accept )
                    kappa(i, j) = newkappa;
                end
                
                logkappaq(i, j) = logkappaq(i, j) + 1 / it^0.55 * (accept - 0.3);
                
                if ( (isnan(kappa(i, j)) == 1) || (isreal(kappa(i, j)) ==0) )
                    stop;
                end
            end
        end
        
        
        for i = 1:T
            chi1 = (target(i) - sum(data(i, :) .* beta(i, :)))^2;
            if ( i == 1 )
                lam1 = kappasigmasq(1) + lambdasigma - 0.5;
                psi1 = 2 * (lambdasigma / musigma + rhosigma / (1 - rhosigma) * lambdasigma / musigma);
            elseif ( i == T )
                lam1 = kappasigmasq(i - 1) + lambdasigma - 0.5;
                psi1 = 2 * (lambdasigma / musigma + rhosigma / (1 - rhosigma) * lambdasigma / musigma);
            else
                lam1 = kappasigmasq(i) + kappasigmasq(i - 1) + lambdasigma - 0.5;
                psi1 = 2 * (lambdasigma / musigma + 2 * rhosigma / (1 - rhosigma) * lambdasigma / musigma);
            end
            
            newsigmasq = sigmasq(i) * exp(exp(logsigmasqsd(i)) * randn);
            
            loglike = (lam1 - 1) * log(sigmasq(i)) - 0.5 * chi1 ./ sigmasq(i) - 0.5 * psi1 * sigmasq(i);
            newloglike = (lam1 - 1) * log(newsigmasq) - 0.5 * chi1 ./ newsigmasq - 0.5 * psi1 * newsigmasq;
            
            logaccept = newloglike - loglike + log(newsigmasq) - log(sigmasq(i));
            
            accept = 1;
            if ( (isnan(logaccept) == 1) || (isinf(logaccept) == 1) )
                accept = 0;
            elseif ( logaccept < 0 )
                accept = exp(logaccept);
            end
            
            sigmasq1accept(i) = sigmasq1accept(i) + accept;
            sigmasq1count(i) = sigmasq1count(i) + 1;
            
            if ( rand < accept )
                sigmasq(i) = newsigmasq;
            end
            
            logsigmasqsd(i) = logsigmasqsd(i) + 1 / it^0.55 * (accept - 0.3);
        end
        
        for i = 1:(T-1)
            newkappasigmasq = kappasigmasq(i) + (2 * ceil(rand) - 1) * geornd(1 / (1 + exp(logkappasigmasqq(i))));

            if ( newkappasigmasq < 0 )
                accept = 0;
            else
                lam1 = lambdasigma + kappasigmasq(i);
                gam1 = lambdasigma / musigma + rhosigma / (1 - rhosigma) * lambdasigma / musigma;
                loglike = lam1 * log(gam1) - gammaln(lam1) + (lam1 - 1) * log(sigmasq(i + 1));
                pnmean = sigmasq(i) * rhosigma / (1 - rhosigma) * lambdasigma / musigma;
                loglike = loglike + kappasigmasq(i) * log(pnmean) - gammaln(kappasigmasq(i) + 1);
                
                lam1 = lambdasigma + newkappasigmasq;
                gam1 = lambdasigma / musigma + rhosigma / (1 - rhosigma) * lambdasigma / musigma;
                newloglike = lam1 * log(gam1) - gammaln(lam1) + (lam1 - 1) * log(sigmasq(i + 1));
                pnmean = sigmasq(i) * rhosigma / (1 - rhosigma) * lambdasigma / musigma;
                newloglike = newloglike + newkappasigmasq * log(pnmean) - gammaln(newkappasigmasq + 1);
                
                logaccept = newloglike - loglike;
                
                accept = 1;
                if ( (isnan(logaccept) == 1) || (isinf(logaccept) == 1) )
                    accept = 0;
                elseif ( logaccept < 0 )
                    accept = exp(logaccept);
                end
            end
            
            kappalambdasigmaccept = kappalambdasigmaccept + accept;
            kappasigmasqcount = kappasigmasqcount + 1;
            
            if ( rand < accept )
                kappasigmasq(i) = newkappasigmasq;
            end
            
            logkappasigmasqq(i) = logkappasigmasqq(i) + 1 / it^0.55 * (accept - 0.3);
            
            
            if ( (isnan(kappasigmasq(i)) == 1) || (isreal(kappasigmasq(i)) ==0) )
                stop;
            end
        end
    end
    
    zstar = rand(1, p) < (5 / p);
    
    
    targetstar = target - (sum(data(:, zstar==0) .* beta(:, zstar==0), 2))';
    datastar = data(:, zstar==1);
    Psistar = Psi(:, zstar==1);
    kappastar = kappa(:, zstar==1);
    
    [meankf, varkf, loglike] = kf_NGAR(datastar, targetstar, Psistar, sigmasq, rhobeta(zstar==1));
    

    
    
    xstar = [log(lambdasigma); log(musigma); log(rhosigma) - log(1 - rhosigma)];
    
    if ( it < 100 )
        newxstar = xstar + [exp(loglambdasigmasd); exp(loggammasigmasqsd); exp(logrhosigmasd)] .* randn(3, 1);
    else
        varstar1 = ([sum2sigmasq(1) sum2sigmasq(2) sum2sigmasq(4); sum2sigmasq(2) sum2sigmasq(3) sum2sigmasq(5); sum2sigmasq(4) sum2sigmasq(5) sum2sigmasq(6)] - sum1sigmasq * sum1sigmasq' / it) / (it - 1);
        
        newxstar = xstar + chol(exp(logscalesigmasq) * varstar1)' * randn(3, 1);
    end
    
    newlambdasigma = exp(newxstar(1));
    newmusigma = exp(newxstar(2));
    newrhosigma = exp(newxstar(3)) / (1 + exp(newxstar(3)));
    
    if ( newrhosigma > limit1 )
        accept = 0;
    else
        newsigmasq = sigmasq;
        newkappasigmasq = kappasigmasq;
        
        newsigmasq(1) = sigmasq(1) * newmusigma / musigma;
        if ( newlambdasigma > lambdasigma )
            newsigmasq(1) = newsigmasq(1) + gamrnd(newlambdasigma - lambdasigma, newmusigma / newlambdasigma);
        else
            newsigmasq(1) = newsigmasq(1) * betarnd(newlambdasigma, lambdasigma - newlambdasigma);
        end
        for i = 2:T
            oldmean = rhosigma / (1 - rhosigma) * lambdasigma / musigma * sigmasq(i - 1);
            newmean = newrhosigma / (1 - newrhosigma) * newlambdasigma / newmusigma * newsigmasq(i - 1);
            
            if (newmean > oldmean)
                newkappasigmasq(i - 1) = kappasigmasq(i - 1) + poissrnd(newmean - oldmean);
            else
                newkappasigmasq(i - 1) = binornd(kappasigmasq(i - 1), newmean / oldmean);
            end
            
            oldlam = kappasigmasq(i - 1) + lambdasigma;
            oldgam = rhosigma / (1 - rhosigma) * lambdasigma / musigma + lambdasigma / musigma;
            newlam = newkappasigmasq(i - 1) + newlambdasigma;
            newgam = newrhosigma / (1 - newrhosigma) * newlambdasigma / newmusigma + newlambdasigma / newmusigma;
            
            newsigmasq(i) = sigmasq(i) * oldgam / newgam;
            
            if (newlam > oldlam)
                newsigmasq(i) = newsigmasq(i) + gamrnd(newlam - oldlam, 1 / newgam);
            else
                newsigmasq(i) = newsigmasq(i) * betarnd(newlam, oldlam - newlam);
            end
        end
        
        [newmeankf, newvarkf, newloglike] = kf_NGAR(datastar, targetstar, Psistar, newsigmasq, rhobeta(zstar==1));
        logaccept = newloglike - loglike + 3 * (log(newlambdasigma) - log(lambdasigma)) - 1 * (newlambdasigma - lambdasigma);
        logaccept = logaccept + log(newmusigma) - log(musigma);
        logaccept = logaccept - (1 + 0.5) * log(1 + newmusigma) + (1 + 0.5) * log(1 + musigma);
        logaccept = logaccept + log(1 / rhosigma + 1 / (1 - rhosigma)) - log(1 / newrhosigma + 1 / (1 - newrhosigma));
        logaccept = logaccept + 40*0.95*(log(newrhosigma) - log(rhosigma)) + 40*0.05*(log(1 - newrhosigma) - log(1 - rhosigma));
        
        accept = 1;
        if ( (isnan(logaccept) == 1) || (isinf(logaccept) == 1) )
            accept = 0;
        elseif ( logaccept < 0 )
            accept = exp(logaccept);
        end
    end
    sigmasqparamaccept = sigmasqparamaccept + accept;
    sigmasqparamcount = sigmasqparamcount + 1;
    
    if ( rand < accept )
        lambdasigma = newlambdasigma;
        musigma = newmusigma;
        rhosigma = newrhosigma;
        sigmasq = newsigmasq;
        kappasigmasq = newkappasigmasq;
        loglike = newloglike;
        meankf = newmeankf;
        varkf = newvarkf;
    end
    
    if ( it < 100 )
        loglambdasigmasd = loglambdasigmasd + 1 / it^0.55 * (accept - 0.3);
        loggammasigmasqsd = loggammasigmasqsd + 1 / it^0.55 * (accept - 0.3);
        logrhosigmasd = logrhosigmasd + 1 / it^0.55 * (accept - 0.3);
    else
        logscalesigmasq = logscalesigmasq + 1 / (it - 99)^0.55 * (accept - 0.3);
    end
    
    x1 = log(lambdasigma);
    x2 = log(musigma);
    x3 = log(rhosigma) - log(1 - rhosigma);
    
    sum1sigmasq(1) = sum1sigmasq(1) + x1;
    sum1sigmasq(2) = sum1sigmasq(2) + x2;
    sum1sigmasq(3) = sum1sigmasq(3) + x3;
    
    sum2sigmasq(1) = sum2sigmasq(1) + x1^2;
    sum2sigmasq(2) = sum2sigmasq(2) + x1*x2;
    sum2sigmasq(3) = sum2sigmasq(3) + x2^2;
    sum2sigmasq(4) = sum2sigmasq(4) + x1*x3;
    sum2sigmasq(5) = sum2sigmasq(5) + x2*x3;
    sum2sigmasq(6) = sum2sigmasq(6) + x3^2;
    
    


    
    
    counter = 0;
    for j = 1:p
        if ( zstar(j) == 1 )
            counter = counter + 1;
            
            xstar = [log(lambda(j)); log(mu(j)); log(rho(j))-log(1-rho(j)); log(rhobeta(j))-log(1-rhobeta(j))];
            
            if ( it < start_adap )
                newxstar = xstar + [exp(loglambdasd(j)); exp(logmeansd(j)); exp(logrhosd(j)); exp(logrhobetasd(j))] .* randn(4, 1);
            else
                sxx = [sum2(1, j) sum2(2, j) sum2(4, j) sum2(7, j); sum2(2, j) sum2(3, j) sum2(5, j) sum2(8, j); sum2(4, j) sum2(5, j) sum2(6, j) sum2(9, j); sum2(7, j) sum2(8, j) sum2(9, j) sum2(10, j)];
                varstar1 = (sxx - sum1(:, j) * sum1(:, j)' / (it - start_samples)) / (it - start_samples - 1);
                
                newxstar = xstar + chol(exp(logscale(j)) * varstar1)' * randn(4, 1);
            end
            
            newlambda = exp(newxstar(1));
            newmu = exp(newxstar(2));
            newrho = exp(newxstar(3)) / (1 + exp(newxstar(3)));
            newrhobeta = rhobeta;
            newrhobeta(j) = exp(newxstar(4)) / (1 + exp(newxstar(4)));
            newDelta = newrho / (1 - newrho) * newlambda / newmu;
            
            if ( (newrhobeta(j) > limit1) || (newrho > limit1) )
                accept = 0;
            else
                newPsistar = Psistar;
                newkappastar = kappastar;
                
                newPsistar(1, counter) = Psistar(1, counter) * newmu / mu(j);
                if ( newlambda > lambda(j) )
                    newPsistar(1, counter) = newPsistar(1, counter) + gamrnd(newlambda - lambda(j), newmu / newlambda);
                else
                    newPsistar(1, counter) = newPsistar(1, counter) * betarnd(newlambda, lambda(j) - newlambda);
                end
                for i = 2:T
                    oldmean = Delta(j) * Psistar(i - 1, counter);
                    newmean = newDelta * newPsistar(i - 1, counter);
                    
                    if (newmean > oldmean)
                        newkappastar(i - 1, counter) = kappastar(i - 1, counter) + poissrnd(newmean - oldmean);
                    else
                        newkappastar(i - 1, counter) = binornd(kappastar(i - 1, counter), newmean / oldmean);
                    end
                    
                    oldlam = kappastar(i - 1, counter) + lambda(j);
                    oldgam = Delta(j) + lambda(j) / mu(j);
                    newlam = newkappastar(i - 1, counter) + newlambda;
                    newgam = newDelta + newlambda / newmu;
                    
                    newPsistar(i, counter) = Psistar(i, counter) * oldgam / newgam;
                    
                    if (newlam > oldlam)
                        newPsistar(i, counter) = newPsistar(i, counter) + gamrnd(newlam - oldlam, 1 / newgam);
                    else
                        newPsistar(i, counter) = newPsistar(i, counter) * betarnd(newlam, oldlam - newlam);
                    end
                end
                
                [newmeankf, newvarkf, newloglike] = kf_NGAR(datastar, targetstar, newPsistar, sigmasq, newrhobeta(zstar==1));
                logaccept = newloglike - loglike;
                if (j == 1)
                    logaccept = logaccept + 2 * log(newlambda) - 2 * log(lambda(j)) - 4 * log(0.5 + newlambda) + 4 * log(0.5 + lambda(j));
                    logaccept = logaccept + log(newmu) - log(mu(j));
                    logaccept = logaccept - (1 + 0.5) * log(1 + newmu) + (1 + 0.5) * log(1 + mu(j));
                else
                    logaccept = logaccept + 2 * log(newlambda) - 2 * log(lambda(j)) - 4 * log(0.5 + newlambda) + 4 * log(0.5 + lambda(j));
                    logaccept = logaccept + lambdastar * (log(newmu) - log(mu(j))) - lambdastar / mustar * (newmu - mu(j));
                end
                logaccept = logaccept + log(1 / rho(j) + 1 / (1 - rho(j))) - log(1 / newrho + 1 / (1 - newrho));
                logaccept = logaccept + 80*0.97*(log(newrho) - log(rho(j))) + 80*0.03*(log(1 - newrho) - log(1 - rho(j)));
                logaccept = logaccept + log(1 / rhobeta(j) + 1 / (1 - rhobeta(j))) - log(1 / newrhobeta(j) + 1 / (1 - newrhobeta(j)));
                logaccept = logaccept + 80*0.97*(log(newrhobeta(j)) - log(rhobeta(j))) + 80*0.03*(log(1 - newrhobeta(j)) - log(1 - rhobeta(j)));
                
                accept = 1;
                if ( (isnan(logaccept) == 1) || (isinf(logaccept) == 1) )
                    accept = 0;
                elseif ( logaccept < 0 )
                    accept = exp(logaccept);
                end
            end
            
            Psiparamaccept(j) = Psiparamaccept(j) + accept;
            Psiparamcount(j) = Psiparamcount(j) + 1;
            
            if ( rand < accept )
                lambda(j) = newlambda;
                mu(j) = newmu;
                rho(j) = newrho;
                rhobeta(j) = newrhobeta(j);
                Delta(j) = newDelta;
                Psistar(:, counter) = newPsistar(:, counter);
                kappastar(:, counter) = newkappastar(:, counter);
                loglike = newloglike;
                meankf = newmeankf;
                varkf = newvarkf;
                
            end
            
            if ( it < start_adap )
                loglambdasd(j) = loglambdasd(j) + 1 / it^0.55 * (accept - 0.3);
                logmeansd(j) = logmeansd(j) + 1 / it^0.55 * (accept - 0.3);
                logrhosd(j) = logrhosd(j) + 1 / it^0.55 * (accept - 0.3);
                logrhobetasd(j) = logrhobetasd(j) + 1 / it^0.55 * (accept - 0.3);
            else
                logscale(j) = logscale(j) + 1 / (it - 99)^0.55 * (accept - 0.3);
            end
        end
    end
    
    if ( it >= start_samples )
        x1 = log(lambda);
        x2 = log(mu);
        x3 = log(rho) - log(1 - rho);
        x4 = log(rhobeta) - log(1 - rhobeta);
        
        sum1(1, :) = sum1(1, :) + x1;
        sum1(2, :) = sum1(2, :) + x2;
        sum1(3, :) = sum1(3, :) + x3;
        sum1(4, :) = sum1(4, :) + x4;
        
        sum2(1, :) = sum2(1, :) + x1.^2;
        sum2(2, :) = sum2(2, :) + x1 .* x2;
        sum2(3, :) = sum2(3, :) + x2.^2;
        sum2(4, :) = sum2(4, :) + x1 .* x3;
        sum2(5, :) = sum2(5, :) + x2 .* x3;
        sum2(6, :) = sum2(6, :) + x3.^2;
        sum2(7, :) = sum2(7, :) + x1 .* x4;
        sum2(8, :) = sum2(8, :) + x2 .* x4;
        sum2(9, :) = sum2(9, :) + x3 .* x4;
        sum2(10, :) = sum2(10, :) + x4.^2;
    end
    
    
    
    Psi(:, zstar == 1) = Psistar;
    kappa(:, zstar == 1) = kappastar;
    
    newbeta = zeros(T, sum(zstar));
    
    checkstar = 1;
    [cholstar, check] = chol(varkf(:, :, T));
    if ( check == 0 )
        newbeta(T, :) = (meankf(:, T) + cholstar' * randn(size(cholstar, 2), 1))';
    else
        checkstar = 0;
    end
    for i = (T-1):-1:1
        Gkal = diag(rhobeta(zstar == 1) .* sqrt(Psistar(i + 1, :) ./ Psistar(i, :)));
        invQ = diag(1 ./ (1 - rhobeta(zstar == 1).^2) .* 1./Psistar(i + 1, :));
        invvarkf = inv(varkf(:, :, i));
        varfb = inv(invvarkf + Gkal' * invQ * Gkal);
        meanfb = varfb * (invvarkf * meankf(:, i) + Gkal' * invQ * newbeta(i + 1, :)');
        [cholstar, check] = chol(varfb);
        if ( check == 0 )
            newbeta(i, :) = (meanfb + cholstar' * randn(size(cholstar, 2), 1))';
        else
            checkstar = 0;
        end
    end
    
    if ( (checkstar == 1) && (sum(sum(isnan(newbeta))) == 0) )
        beta(:, zstar==1) = newbeta;
    end
    
    
    
    
    
    % Update mustar
    newmustar = mustar * exp(mugammasd * randn);
       
    logaccept = (p - 1) * lambdastar * (log(mustar) - log(newmustar));
    logaccept = logaccept - lambdastar * (1 / newmustar - 1 / mustar) * sum(mu(2:end));
    logaccept = logaccept + log(newmustar) - log(mustar) - 3 * log(newmustar + mumean1b) + 3 * log(mustar + mumean1b);
    
    accept = 1;
    if ( (isnan(logaccept) == 1) || (isinf(logaccept) == 1) )
        accept = 0;
    elseif ( logaccept < 0 )
        accept = exp(logaccept);
    end
    
    
    mugammaaccept = mugammaaccept + accept;
    mugammacount = mugammacount + 1;
    
    if ( rand < accept )
        mustar = newmustar;
    end
    
    newmugammasd = mugammasd + 1 / it^0.5 * (accept - 0.3);
    
    if ( (newmugammasd > 10^(-3)) && (newmugammasd < 10^3) )
        mugammasd = newmugammasd;
    end
    
    
    
    
    % Update lambdastar
    newlambdastar = lambdastar * exp(vgammasd * randn);
    
    logaccept = (p - 1) * (newlambdastar * log(newlambdastar / mustar) - lambdastar * log(lambdastar / mustar));
    logaccept = logaccept - (p - 1) * (gammaln(newlambdastar) - gammaln(lambdastar));
    logaccept = logaccept + (newlambdastar - lambdastar) * sum(log(mu(2:end)));
    logaccept = logaccept - (newlambdastar - lambdastar) / mustar * sum(mu(2:end));
    logaccept = logaccept + log(newlambdastar) - log(lambdastar) - 1 / mulambdastar * (newlambdastar - lambdastar);
    
    accept = 1;
    if ( (isnan(logaccept) == 1) || (isinf(logaccept) == 1) )
        accept = 0;
    elseif ( logaccept < 0 )
        accept = exp(logaccept);
    end
    
    vgammaaccept = vgammaaccept + accept;
    vgammacount = vgammacount + 1;
    
    if ( rand < accept )
        lambdastar = newlambdastar;
    end
    
    newvgammasd = vgammasd + 1 / it^0.5 * (accept - 0.3);
    
    if ( (newvgammasd > 10^(-3)) && (newvgammasd < 10^3) )
        vgammasd = newvgammasd;
    end
    
    
    if ( (it > burnin) && (mod(it - burnin, every) == 0) )
        holdbeta(:, :, (it - burnin) / every) = beta;
        holdPsi(:, :, (it - burnin) / every) = Psi;
        holdsigmasq(:, (it - burnin) / every) = sigmasq;
        holdlambda(:, (it - burnin) / every) = lambda;
        holdmu(:, (it - burnin) / every) = mu;
        holdrho(:, (it - burnin) / every) = rho;
        holdrhobeta(:, (it - burnin) / every) = rhobeta;
        holdlambdasigma((it - burnin) / every) = lambdasigma;
        holdmusigma((it - burnin) / every) = musigma;
        holdrhosigma((it - burnin) / every) = rhosigma;
        holdlambdastar((it - burnin) / every) = lambdastar;
        holdmustar((it - burnin) / every) = mustar;
    end
end

beta = holdbeta;
sigma = holdsigmasq;

fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8)
% fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8)
toc; % Stop timer and print total time

% output = struct('beta', holdbeta, 'sigmasq', holdsigmasq, 'Psi', holdPsi, 'lambda', holdlambda, 'mu', holdmu, 'rhobeta', holdrhobeta, 'rho', holdrho, 'lambdasigma', holdlambdasigma, 'musigma', holdmusigma, 'rhosigma', holdrhosigma, 'lambdastar', holdlambdastar, 'mustar', holdmustar);