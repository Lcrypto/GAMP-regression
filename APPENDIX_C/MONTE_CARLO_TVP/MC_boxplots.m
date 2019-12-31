figure(1)
load('MC1_T30.mat')
MAD_GAMP = mean(abs(BETA(:,:,1,2) - BETA(:,:,1,1)),2);
MAD_TVP1 =  mean(abs(BETA(:,:,1,3) - BETA(:,:,1,1)),2);
MAD_TVP2 =  mean(abs(BETA(:,:,1,4) - BETA(:,:,1,1)),2);

MAD = [MAD_GAMP, MAD_TVP1, MAD_TVP2];
subplot(1,3,1)
boxplot(MAD)
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'TVP-GAMP', 'MCMC (tight-prior)', 'MCMC (loose-prior)'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('T=30')

load('MC1_T200.mat')
MAD_GAMP = mean(abs(BETA(:,:,1,2) - BETA(:,:,1,1)),2);
MAD_TVP1 =  mean(abs(BETA(:,:,1,3) - BETA(:,:,1,1)),2);
MAD_TVP2 =  mean(abs(BETA(:,:,1,4) - BETA(:,:,1,1)),2);

MAD = [MAD_GAMP, MAD_TVP1, MAD_TVP2];
subplot(1,3,2)
boxplot(MAD)%,'orientation','horizontal')
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'TVP-GAMP', 'MCMC (tight-prior)', 'MCMC (loose-prior)'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('T=200')

load('MC1_T500.mat')
MAD_GAMP = mean(abs(BETA(:,:,1,2) - BETA(:,:,1,1)),2);
MAD_TVP1 =  mean(abs(BETA(:,:,1,3) - BETA(:,:,1,1)),2);
MAD_TVP2 =  mean(abs(BETA(:,:,1,4) - BETA(:,:,1,1)),2);

MAD = [MAD_GAMP, MAD_TVP1, MAD_TVP2];
subplot(1,3,3)
boxplot(MAD)%,'orientation','horizontal')
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'TVP-GAMP', 'MCMC (tight-prior)', 'MCMC (loose-prior)'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('T=500')


figure(2)
load('MC1alt_T30.mat')
MAD_GAMP = mean(abs(BETA(:,:,1,2) - BETA(:,:,1,1)),2);
MAD_TVP1 =  mean(abs(BETA(:,:,1,3) - BETA(:,:,1,1)),2);
MAD_TVP2 =  mean(abs(BETA(:,:,1,4) - BETA(:,:,1,1)),2);

MAD = [MAD_GAMP, MAD_TVP1, MAD_TVP2];
subplot(1,3,1)
boxplot(MAD)
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'TVP-GAMP', 'MCMC (tight-prior)', 'MCMC (loose-prior)'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('T=30')

load('MC1alt_T200.mat')
MAD_GAMP = mean(abs(BETA(:,:,1,2) - BETA(:,:,1,1)),2);
MAD_TVP1 =  mean(abs(BETA(:,:,1,3) - BETA(:,:,1,1)),2);
MAD_TVP2 =  mean(abs(BETA(:,:,1,4) - BETA(:,:,1,1)),2);

MAD = [MAD_GAMP, MAD_TVP1, MAD_TVP2];
subplot(1,3,2)
boxplot(MAD)%,'orientation','horizontal')
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'TVP-GAMP', 'MCMC (tight-prior)', 'MCMC (loose-prior)'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('T=200')

load('MC1alt_T500.mat')
MAD_GAMP = mean(abs(BETA(:,:,1,2) - BETA(:,:,1,1)),2);
MAD_TVP1 =  mean(abs(BETA(:,:,1,3) - BETA(:,:,1,1)),2);
MAD_TVP2 =  mean(abs(BETA(:,:,1,4) - BETA(:,:,1,1)),2);

MAD = [MAD_GAMP, MAD_TVP1, MAD_TVP2];
subplot(1,3,3)
boxplot(MAD)%,'orientation','horizontal')
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'TVP-GAMP', 'MCMC (tight-prior)', 'MCMC (loose-prior)'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('T=500')


figure(3)
load('MC1alt2_T30.mat')
MAD_GAMP = mean(abs(BETA(:,:,1,2) - BETA(:,:,1,1)),2);
MAD_TVP1 =  mean(abs(BETA(:,:,1,3) - BETA(:,:,1,1)),2);
MAD_TVP2 =  mean(abs(BETA(:,:,1,4) - BETA(:,:,1,1)),2);

MAD = [MAD_GAMP, MAD_TVP1, MAD_TVP2];
subplot(1,3,1)
boxplot(MAD)
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'TVP-GAMP', 'MCMC (tight-prior)', 'MCMC (loose-prior)'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('T=30')

load('MC1alt2_T200.mat')
MAD_GAMP = mean(abs(BETA(:,:,1,2) - BETA(:,:,1,1)),2);
MAD_TVP1 =  mean(abs(BETA(:,:,1,3) - BETA(:,:,1,1)),2);
MAD_TVP2 =  mean(abs(BETA(:,:,1,4) - BETA(:,:,1,1)),2);

MAD = [MAD_GAMP, MAD_TVP1, MAD_TVP2];
subplot(1,3,2)
boxplot(MAD)%,'orientation','horizontal')
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'TVP-GAMP', 'MCMC (tight-prior)', 'MCMC (loose-prior)'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('T=200')

load('MC1alt2_T500.mat')
MAD_GAMP = mean(abs(BETA(:,:,1,2) - BETA(:,:,1,1)),2);
MAD_TVP1 =  mean(abs(BETA(:,:,1,3) - BETA(:,:,1,1)),2);
MAD_TVP2 =  mean(abs(BETA(:,:,1,4) - BETA(:,:,1,1)),2);

MAD = [MAD_GAMP, MAD_TVP1, MAD_TVP2];
subplot(1,3,3)
boxplot(MAD)%,'orientation','horizontal')
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'TVP-GAMP', 'MCMC (tight-prior)', 'MCMC (loose-prior)'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('T=500')
