function [y,x,beta] = ar_dgp(T)

y0 = rand;
% Generate from AR(4) model with similar parameters to those from an
% estimated AR(4) model for US real GDP (1960Q1 to 2016Q4)
beta = [0.4, 0.22, 0.05, 0.14];
y = zeros(T+100,1);
for t = 1:T+100
   if t == 1
       y(t,:) = beta(1)*y0 + randn;       
   elseif t == 2
       y(t,:) = beta(1)*y(t-1,:) + beta(2)*y0 + randn;
   elseif t == 3  
       y(t,:) = beta(1)*y(t-1,:) + beta(2)*y(t-2,:) + beta(3)*y0 + randn;
   elseif t == 4
       y(t,:) = beta(1)*y(t-1,:) + beta(2)*y(t-2,:) + beta(3)*y(t-3,:) + beta(4)*y0 + randn;
   elseif t > 4
       y(t,:) = beta(1)*y(t-1,:) + beta(2)*y(t-2,:) + beta(3)*y(t-3,:) + beta(4)*y(t-4,:) + randn;
   end    
end

x = [y(end-T:end-1,:), y(end-T-1:end-2,:), y(end-T-2:end-3,:), y(end-T-3:end-4,:)];
y = y(end-T+1:end,:);
