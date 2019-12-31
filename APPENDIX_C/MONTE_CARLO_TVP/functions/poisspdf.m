function y = poisspdf(x,lambda)

if x == 0 || x == 1
    y=(lambda^x)*exp(-lambda);   
else
    y=((lambda^x)/factorial(x))*exp(-lambda);
end