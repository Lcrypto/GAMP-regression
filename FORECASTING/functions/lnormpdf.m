
function y = lnormpdf(x,mu,sigma)

if nargin<1
    error(message('stats:normpdf:TooFewInputs'));
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end


% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    y = log(1./(sqrt(2*pi) .* sigma)) + (-0.5 * ((x - mu)./sigma).^2);
catch
    error(message('stats:normpdf:InputSizeMismatch'));
end