% Forecasting using Bayesian Hierarchical Shrinkage regressions
%==========================================================================
% This code replicates the forecasting results in Korobilis (2010)
%==========================================================================
 
% Reset everything, clear memory and screen
clear; close all; clc;

% Add path of random number generators
addpath('functions')
addpath('data')

% ===========================| LOAD DATA |=================================
[A,B,~] = xlsread('FREDMD.xlsx','Data');
data    = A(2:end,:);
tcode   = A(1,:);      % Vector of transformations
Ynames  = B(1,2:end);  % The mnemonic for the variables in the panel
%---------------

tcode(tcode==6)=5;

Xraw = 0*data;
for kk = 1:size(data,2)
    if tcode(kk) == 5 || tcode(kk) == 6
        Xraw(:,kk) = 100*transx(data(:,kk),tcode(kk));
    else
        Xraw(:,kk) = transx(data(:,kk),tcode(kk));
    end
end


focus_variables  = [1,3,6,7,10,19,20,23,24,64,70,87,93,94,95,105];
Xraw = Xraw(:,focus_variables);   Ynames = Ynames(focus_variables);
[T,p] = size(Xraw);

VOL = zeros(p,T,2);
for i = 1:p
    y = Xraw(:,i);
    
    % GAMP estimates
    [beta_SBLTV,~,sigma_SBLTV,~]     = AMPSBL_SV([ones(T,1), eye(T)],y,1000,1e-10,1e-10);

    % GARCH(1,1) estimates
    ToEstMdl = garch(1,1);
    EstMdl = estimate(ToEstMdl,(y - mean(y)));
    R      = infer(EstMdl,y);
    
    % Save estimates in VOL matrix
    VOL(i,:,1) = sigma_SBLTV; VOL(i,:,2) = R;
end
clc;
pause(1);

% Plot estimates for all focus variables
figure
for j = 1:16
    subplot(4,4,j)    
    plot(VOL(j,:,1),'-','LineWidth',2)
    hold all
    plot(VOL(j,:,2),'--','LineWidth',2)
    if j == 1
        legend({'GAMP-SV';'GARCH(1,1)'},'Location','northwest')
    end
    title([cell2mat(Ynames(j))])
end

% 
% %% Prior sensitivity code, estimate a trend for inflation (CPI, Total, All items)
% y = log(data(13:end,95)./data(1:end-12,95)); % This is CPIAUSL
% T = size(y,1)-1;
% 
% figure
% subplot(3,1,1)
% [beta_SBLTV,~,~,~]     = AMPSBL_SV([ones(T,1), eye(T), y(1:end-1)],y(2:end),1000,1,1e-10);
% plot(y,'-','LineWidth',2);
% hold all;
% plot((beta_SBLTV(1)+beta_SBLTV(2:end-1)),'--','LineWidth',2);
% legend({'\pi_{t}';'\tau_{t}'},'Location','northwest')
% title('(a)')
% 
% subplot(3,1,2)
% [beta_SBLTV,~,~,~]     = AMPSBL_SV([ones(T,1), eye(T), y(1:end-1)],y(2:end),1000,1,0.1);
% plot(y,'-','LineWidth',2);
% hold all;
% plot((beta_SBLTV(1)+beta_SBLTV(2:end-1)),'--','LineWidth',2);
% legend({'\pi_{t}';'\tau_{t}'},'Location','northwest')
% title('(b)')
% 
% subplot(3,1,3)
% [beta_SBLTV,~,~,~]     = AMPSBL_SV([ones(T,1), eye(T), y(1:end-1)],y(2:end),1000,1,10);
% plot(y,'-','LineWidth',2);
% hold all;
% plot((beta_SBLTV(1)+beta_SBLTV(2:end-1)),'--','LineWidth',2);
% legend({'\pi_{t}';'\tau_{t}'},'Location','northwest')
% title('(c)')
% 
% 
% 
