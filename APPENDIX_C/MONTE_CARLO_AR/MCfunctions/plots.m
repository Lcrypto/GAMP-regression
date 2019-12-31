% Plot some coefficients
figure
subplot(2,2,1)
plot(beta_SBLTV(1:p+q:24750),'LineWidth',2)
hold all
plot(squeeze(mean(beta_TVPAR(:,1,:),3)),'LineWidth',2)
legend({'GAMP-SNS','TVP-AR'})
subplot(2,2,2)
plot(beta_SBLTV(2:p+q:24750),'LineWidth',2)
hold all
plot(squeeze(mean(beta_TVPAR(:,2,:),3)),'LineWidth',2)
legend({'GAMP-SNS','TVP-AR'})
subplot(2,2,3)
plot(beta_SBLTV(3:p+q:24750),'LineWidth',2)
hold all
plot(squeeze(mean(beta_TVPAR(:,3,:),3)),'LineWidth',2)
legend({'GAMP-SNS','TVP-AR'})


% Plot some coefficients
figure
subplot(2,2,1)
plot(beta_SSVSTV(1:q:330),'LineWidth',2)
hold all
plot(squeeze(mean(beta_TVPAR(:,1,:),3)),'LineWidth',2)
legend({'GAMP-SNS','TVP-AR'})
subplot(2,2,2)
plot(beta_SSVSTV(2:q:330),'LineWidth',2)
hold all
plot(squeeze(mean(beta_TVPAR(:,2,:),3)),'LineWidth',2)
legend({'GAMP-SNS','TVP-AR'})
subplot(2,2,3)
plot(beta_SSVSTV(3:q:330),'LineWidth',2)
hold all
plot(squeeze(mean(beta_TVPAR(:,3,:),3)),'LineWidth',2)
legend({'GAMP-SNS','TVP-AR'})


% Plot some coefficients
figure
subplot(2,2,1)
plot(beta_SBLTV(1:q:330),'LineWidth',2)
hold all
plot(squeeze(mean(beta_TVPAR(:,1,:),3)),'LineWidth',2)
legend({'GAMP-SNS','TVP-AR'})
subplot(2,2,2)
plot(beta_SBLTV(2:q:330),'LineWidth',2)
hold all
plot(squeeze(mean(beta_TVPAR(:,2,:),3)),'LineWidth',2)
legend({'GAMP-SNS','TVP-AR'})
subplot(2,2,3)
plot(beta_SBLTV(3:q:330),'LineWidth',2)
hold all
plot(squeeze(mean(beta_TVPAR(:,3,:),3)),'LineWidth',2)
legend({'GAMP-SNS','TVP-AR'})