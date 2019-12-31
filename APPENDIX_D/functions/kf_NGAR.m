function [meankf, varkf, loglike] = kf_NGAR(data, y, Psi, sigmasq, rhobeta)

p = size(data, 2);
T = length(y);

meankf=zeros(p, T);
varkf=zeros(p, p, T);

aminus = zeros(p, 1);
Pminus = diag(Psi(1, :));
datastar = data(1, :);
e = y(1) - datastar * aminus;
invF = inv(sigmasq(1) + datastar * Pminus * datastar');
meankf(:, 1) = aminus + Pminus * datastar' * invF * e;
varkf(:, :, 1) = Pminus - Pminus * datastar' * invF * datastar * Pminus;
loglike = - 0.5 * e^2 * invF + 0.5 * log(invF);

for i = 2:T
    datastar = data(i, :);
    Q = diag((1 - rhobeta.^2) .* Psi(i, :));
    Gkal = diag(rhobeta .* sqrt(Psi(i, :) ./ Psi(i - 1, :)));
    aminus = Gkal * meankf(:, i - 1);
    Pminus = Gkal * varkf(:, :, i - 1) * Gkal' + Q;
    e = y(i) - datastar * aminus;
    invF = 1 / (sigmasq(i) + (datastar * Pminus )* datastar');    
    PinvF = (Pminus * datastar') * invF;
    meankf(:, i) = aminus + PinvF * e;
    varkf(:, :, i) = Pminus - (PinvF * datastar) * Pminus;
    loglike = loglike - 0.5 * e^2 * invF + 0.5 * log(invF);
end
