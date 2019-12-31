figure
load('MC_AR_30.mat')
mad_GAMP = mean(abs(BETA(:,:,2) - BETA(:,:,1)),2);
mad_OLS = mean(abs(BETA(:,:,3) - BETA(:,:,1)),2);

subplot(1,3,1)
boxplot([mad_GAMP,mad_OLS])
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'GAMP', 'OLS'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('T=30')


load('MC_AR_100.mat')
mad_GAMP = mean(abs(BETA(:,:,2) - BETA(:,:,1)),2);
mad_OLS = mean(abs(BETA(:,:,3) - BETA(:,:,1)),2);

subplot(1,3,2)
boxplot([mad_GAMP,mad_OLS])
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'GAMP', 'OLS'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('T=100')


load('MC_AR_500.mat')
mad_GAMP = mean(abs(BETA(:,:,2) - BETA(:,:,1)),2);
mad_OLS = mean(abs(BETA(:,:,3) - BETA(:,:,1)),2);

subplot(1,3,3)
boxplot([mad_GAMP,mad_OLS])
grid on
set(gca,'fontsize',14)
set(gca,'XTickLabel',{'GAMP', 'OLS'})
fix_xticklabels(gca,0.1,{'FontSize',14})
title('T=500')