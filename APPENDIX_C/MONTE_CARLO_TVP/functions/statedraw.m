% Procedure to draw the indexes of 7-component mixture approximation to the chi-2 density  

for i = 1:t       
    for j = 1:numel(m_s)
        temp1 = (1/sqrt(2*pi*u2_s(j)))*exp(-.5*(((yss(1,i) - Sigtdraw(i) - m_s(j) + 1.2704)^2)/u2_s(j)));
        prwS(j,1) = q_s(j,1)*temp1;
    end
    prwS = prwS./sum(prwS);
    cprwS = cumsum(prwS);
    trand = rand(1,1);
    if trand < cprwS(1,1); imix=1;
    elseif trand < cprwS(2,1), imixS=2;
    elseif trand < cprwS(3,1), imixS=3;
    elseif trand < cprwS(4,1), imixS=4;
    elseif trand < cprwS(5,1), imixS=5;
    elseif trand < cprwS(6,1), imixS=6;
    else imixS=7; 
    end
    statedrawS(i,1)=imixS;
end