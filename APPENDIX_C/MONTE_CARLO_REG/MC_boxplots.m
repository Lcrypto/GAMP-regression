figure(1)
load('MONTE_CARLO_50_100_0.3.mat')
MAD_SBL = mean(abs(BETA(:,:,3)  - BETA(:,:,1)),2);
MAD_LASSO = mean(abs(BETA(:,:,5)  - BETA(:,:,1)),2);
MAD_SNSg = mean(abs(BETA(:,:,6)  - BETA(:,:,1)),2);

MAD = [MAD_SBL, MAD_LASSO, MAD_SNSg];
subplot(1,3,1)
boxplot(MAD)
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'GAMP (SBL)', 'MCMC (LASSO)', 'MCMC (SSVS)'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('p=100')

load('MONTE_CARLO_50_200_0.3.mat')
MAD_SBL = mean(abs(BETA(:,:,3)  - BETA(:,:,1)),2);
MAD_LASSO = mean(abs(BETA(:,:,5)  - BETA(:,:,1)),2);
MAD_SNSg = mean(abs(BETA(:,:,6)  - BETA(:,:,1)),2);

MAD = [MAD_SBL, MAD_LASSO, MAD_SNSg];
subplot(1,3,2)
boxplot(MAD)%,'orientation','horizontal')
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'GAMP (SBL)', 'MCMC (LASSO)', 'MCMC (SSVS)'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('p=200')

load('MONTE_CARLO_50_500_0.3.mat')
MAD_SBL = mean(abs(BETA(:,:,3)  - BETA(:,:,1)),2);
MAD_LASSO = mean(abs(BETA(:,:,5)  - BETA(:,:,1)),2);
MAD_SNSg = mean(abs(BETA(:,:,6)  - BETA(:,:,1)),2);

MAD = [MAD_SBL, MAD_LASSO, MAD_SNSg];
subplot(1,3,3)
boxplot(MAD)%,'orientation','horizontal')
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'GAMP (SBL)', 'MCMC (LASSO)', 'MCMC (SSVS)'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('p=500')


figure(2)
load('MONTE_CARLO_200_100_0.01.mat')
MAD_SBL = mean(abs(BETA(:,:,3)  - BETA(:,:,1)),2);
MAD_LASSO = mean(abs(BETA(:,:,5)  - BETA(:,:,1)),2);
MAD_SNSg = mean(abs(BETA(:,:,6)  - BETA(:,:,1)),2);

MAD = [MAD_SBL, MAD_LASSO, MAD_SNSg];
subplot(1,3,1)
boxplot(MAD)
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'GAMP (SBL)', 'MCMC (LASSO)', 'MCMC (SSVS)'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('p=100')

load('MONTE_CARLO_200_200_0.01.mat')
MAD_SBL = mean(abs(BETA(:,:,3)  - BETA(:,:,1)),2);
MAD_LASSO = mean(abs(BETA(:,:,5)  - BETA(:,:,1)),2);
MAD_SNSg = mean(abs(BETA(:,:,6)  - BETA(:,:,1)),2);

MAD = [MAD_SBL, MAD_LASSO, MAD_SNSg];
subplot(1,3,2)
boxplot(MAD)
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'GAMP (SBL)', 'MCMC (LASSO)', 'MCMC (SSVS)'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('p=200')

load('MONTE_CARLO_200_500_0.01.mat')
MAD_SBL = mean(abs(BETA(:,:,3)  - BETA(:,:,1)),2);
MAD_LASSO = mean(abs(BETA(:,:,5)  - BETA(:,:,1)),2);
MAD_SNSg = mean(abs(BETA(:,:,6)  - BETA(:,:,1)),2);

MAD = [MAD_SBL, MAD_LASSO, MAD_SNSg];
subplot(1,3,3)
boxplot(MAD)
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'GAMP (SBL)', 'MCMC (LASSO)', 'MCMC (SSVS)'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('p=500')


figure(3)
load('MONTE_CARLO_200_100_0.9.mat')
MAD_SBL = mean(abs(BETA(:,:,3)  - BETA(:,:,1)),2);
MAD_LASSO = mean(abs(BETA(:,:,5)  - BETA(:,:,1)),2);
MAD_SNSg = mean(abs(BETA(:,:,6)  - BETA(:,:,1)),2);

MAD = [MAD_SBL, MAD_LASSO, MAD_SNSg];
boxplot(MAD)
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'GAMP (SBL)', 'MCMC (LASSO)', 'MCMC (SSVS)'})
fix_xticklabels(gca,0.1,{'FontSize',14})


% times = [ t_sbl, t_sns, t_lasso, t_snsg]