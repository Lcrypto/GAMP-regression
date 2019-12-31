function i = lookupKt(Kt)

p = length(Kt);
sumK = sum(Kt);
if sumK == 0
    i=1;
elseif sumK == p
    i=p+2;
else
    i = find(Kt==1)  + 1;
end