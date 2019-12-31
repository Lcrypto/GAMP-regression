function [T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13] = create_tables_rolling(sel_series)

if nargin == 0
    sel_series = 1:222;
end

% MSFEs
% Table 1
load('E:\Dropbox\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\rolling\FORECASTING_(_h_=_1_).mat')
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,9,sel_series))),1,9));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MRMSFE([1:8],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); LASSO = entries(:,3); SSVS = entries(:,4); BMA = entries(:,5); BAG = entries(:,6);DFM5 = entries(:,8); OLS = entries(:,7);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T1 = table(SBL,SNS,LASSO,SSVS,BMA,BAG,DFM5,OLS,'RowNames',Name);

% Table 2
load('E:\Dropbox\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\rolling\FORECASTING_(_h_=_2_).mat')
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,9,sel_series))),1,9));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MRMSFE([1:8],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); LASSO = entries(:,3); SSVS = entries(:,4); BMA = entries(:,5); BAG = entries(:,6);DFM5 = entries(:,8); OLS = entries(:,7);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T2 = table(SBL,SNS,LASSO,SSVS,BMA,BAG,DFM5,OLS,'RowNames',Name);

% Table 3
load('E:\Dropbox\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\rolling\FORECASTING_(_h_=_4_).mat')
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,9,sel_series))),1,9));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MRMSFE([1:8],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); LASSO = entries(:,3); SSVS = entries(:,4); BMA = entries(:,5); BAG = entries(:,6);DFM5 = entries(:,8); OLS = entries(:,7);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T3 = table(SBL,SNS,LASSO,SSVS,BMA,BAG,DFM5,OLS,'RowNames',Name);

% Table 4
load('E:\Dropbox\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\rolling\FORECASTING_(_h_=_8_).mat')
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,9,sel_series))),1,9));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MRMSFE([1:8],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); LASSO = entries(:,3); SSVS = entries(:,4); BMA = entries(:,5); BAG = entries(:,6);DFM5 = entries(:,8); OLS = entries(:,7);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T4 = table(SBL,SNS,LASSO,SSVS,BMA,BAG,DFM5,OLS,'RowNames',Name);

% MAFEs
% Table 5
load('E:\Dropbox\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\rolling\FORECASTING_(_h_=_1_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,9,sel_series)),1,9));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MMAFE([1:8],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); LASSO = entries(:,3); SSVS = entries(:,4); BMA = entries(:,5); BAG = entries(:,6);DFM5 = entries(:,8); OLS = entries(:,7);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T5 = table(SBL,SNS,LASSO,SSVS,BMA,BAG,DFM5,OLS,'RowNames',Name);

% Table 6
load('E:\Dropbox\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\rolling\FORECASTING_(_h_=_2_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,9,sel_series)),1,9));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MMAFE([1:8],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); LASSO = entries(:,3); SSVS = entries(:,4); BMA = entries(:,5); BAG = entries(:,6);DFM5 = entries(:,8); OLS = entries(:,7);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T6 = table(SBL,SNS,LASSO,SSVS,BMA,BAG,DFM5,OLS,'RowNames',Name);

% Table 7
load('E:\Dropbox\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\rolling\FORECASTING_(_h_=_4_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,9,sel_series)),1,9));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MMAFE([1:8],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); LASSO = entries(:,3); SSVS = entries(:,4); BMA = entries(:,5); BAG = entries(:,6);DFM5 = entries(:,8); OLS = entries(:,7);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T7 = table(SBL,SNS,LASSO,SSVS,BMA,BAG,DFM5,OLS,'RowNames',Name);

% Table 8
load('E:\Dropbox\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\rolling\FORECASTING_(_h_=_8_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,9,sel_series)),1,9));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MMAFE([1:8],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); LASSO = entries(:,3); SSVS = entries(:,4); BMA = entries(:,5); BAG = entries(:,6);DFM5 = entries(:,8); OLS = entries(:,7);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T8 = table(SBL,SNS,LASSO,SSVS,BMA,BAG,DFM5,OLS,'RowNames',Name);

% Correlations
% Table 9
load('E:\Dropbox\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\rolling\FORECASTING_(_h_=_1_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,9,sel_series)),1,9));
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,9,sel_series))),1,9));  % Mean Relative MSFE (benchmark is ARp)
RMSFE = zeros(8,8);
for i=1:7
    for j=i+1:8
        RMSFE(i,j) = mean(abs(MRMSFE(i,:)-MRMSFE(j,:)));
    end
end
entries = tril(corrcoef(MMAFE(1:8,:)')) + RMSFE - eye(8);
SBL = entries(:,1); SNS = entries(:,2); LASSO = entries(:,3); SSVS = entries(:,4); BMA = entries(:,5); BAG = entries(:,6);DFM5 = entries(:,8); OLS = entries(:,7);
Name = {'SBL','SNS','LASSO','SSVS','BMA','BAG','DFM5','OLS'};
T9 = table(SBL,SNS,LASSO,SSVS,BMA,BAG,DFM5,OLS,'RowNames',Name);

% Hit Rates
% Table 10
load('E:\Dropbox\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\rolling\FORECASTING_(_h_=_1_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,9,sel_series)),1,9));  % Mean Relative MSFE (benchmark is ARp)
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,9,sel_series))),1,9));  % Mean Relative MSFE (benchmark is ARp)
[~,b1] = min(MMAFE(1:8,:)); [~,b2] = min(MRMSFE(1:8,:));
SBL = 100*[sum(b1==1)/length(sel_series);sum(b2==1)/length(sel_series)] ; SNS = 100*[sum(b1==2)/length(sel_series);sum(b2==2)/length(sel_series)]; LASSO = 100*[sum(b1==3)/length(sel_series);sum(b2==3)/length(sel_series)]; SSVS = 100*[sum(b1==4)/length(sel_series);sum(b2==4)/length(sel_series)]; BMA = 100*[sum(b1==5)/length(sel_series);sum(b2==5)/length(sel_series)]; BAG = 100*[sum(b1==6)/length(sel_series);sum(b2==6)/length(sel_series)]; OLS = 100*[sum(b1==7)/length(sel_series);sum(b2==7)/length(sel_series)];DFM5 = 100*[sum(b1==8)/length(sel_series);sum(b2==8)/length(sel_series)];
Name = {'% of lowest MAFE','% of lowest MSFE'};
T10 = table(SBL,SNS,LASSO,SSVS,BMA,BAG,DFM5,OLS,'RowNames',Name);

% Table 11
load('E:\Dropbox\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\rolling\FORECASTING_(_h_=_2_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,9,sel_series)),1,9));  % Mean Relative MSFE (benchmark is ARp)
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,9,sel_series))),1,9));  % Mean Relative MSFE (benchmark is ARp)
[~,b1] = min(MMAFE(1:8,:)); [~,b2] = min(MRMSFE(1:8,:));
SBL = 100*[sum(b1==1)/length(sel_series);sum(b2==1)/length(sel_series)] ; SNS = 100*[sum(b1==2)/length(sel_series);sum(b2==2)/length(sel_series)]; LASSO = 100*[sum(b1==3)/length(sel_series);sum(b2==3)/length(sel_series)]; SSVS = 100*[sum(b1==4)/length(sel_series);sum(b2==4)/length(sel_series)]; BMA = 100*[sum(b1==5)/length(sel_series);sum(b2==5)/length(sel_series)]; BAG = 100*[sum(b1==6)/length(sel_series);sum(b2==6)/length(sel_series)]; OLS = 100*[sum(b1==7)/length(sel_series);sum(b2==7)/length(sel_series)];DFM5 = 100*[sum(b1==8)/length(sel_series);sum(b2==8)/length(sel_series)];
Name = {'% of lowest MAFE','% of lowest MSFE'};
T11 = table(SBL,SNS,LASSO,SSVS,BMA,BAG,DFM5,OLS,'RowNames',Name);

% Table 12
load('E:\Dropbox\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\rolling\FORECASTING_(_h_=_4_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,9,sel_series)),1,9));  % Mean Relative MSFE (benchmark is ARp)
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,9,sel_series))),1,9));  % Mean Relative MSFE (benchmark is ARp)
[~,b1] = min(MMAFE(1:8,:)); [~,b2] = min(MRMSFE(1:8,:));
SBL = 100*[sum(b1==1)/length(sel_series);sum(b2==1)/length(sel_series)] ; SNS = 100*[sum(b1==2)/length(sel_series);sum(b2==2)/length(sel_series)]; LASSO = 100*[sum(b1==3)/length(sel_series);sum(b2==3)/length(sel_series)]; SSVS = 100*[sum(b1==4)/length(sel_series);sum(b2==4)/length(sel_series)]; BMA = 100*[sum(b1==5)/length(sel_series);sum(b2==5)/length(sel_series)]; BAG = 100*[sum(b1==6)/length(sel_series);sum(b2==6)/length(sel_series)]; OLS = 100*[sum(b1==7)/length(sel_series);sum(b2==7)/length(sel_series)];DFM5 = 100*[sum(b1==8)/length(sel_series);sum(b2==8)/length(sel_series)];
Name = {'% of lowest MAFE','% of lowest MSFE'};
T12 = table(SBL,SNS,LASSO,SSVS,BMA,BAG,DFM5,OLS,'RowNames',Name);

% Table 13
load('E:\Dropbox\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\rolling\FORECASTING_(_h_=_8_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,9,sel_series)),1,9));  % Mean Relative MSFE (benchmark is ARp)
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,9,sel_series))),1,9));  % Mean Relative MSFE (benchmark is ARp)
[~,b1] = min(MMAFE(1:8,:)); [~,b2] = min(MRMSFE(1:8,:));
SBL = 100*[sum(b1==1)/length(sel_series);sum(b2==1)/length(sel_series)] ; SNS = 100*[sum(b1==2)/length(sel_series);sum(b2==2)/length(sel_series)]; LASSO = 100*[sum(b1==3)/length(sel_series);sum(b2==3)/length(sel_series)]; SSVS = 100*[sum(b1==4)/length(sel_series);sum(b2==4)/length(sel_series)]; BMA = 100*[sum(b1==5)/length(sel_series);sum(b2==5)/length(sel_series)]; BAG = 100*[sum(b1==6)/length(sel_series);sum(b2==6)/length(sel_series)]; OLS = 100*[sum(b1==7)/length(sel_series);sum(b2==7)/length(sel_series)];DFM5 = 100*[sum(b1==8)/length(sel_series);sum(b2==8)/length(sel_series)];
Name = {'% of lowest MAFE','% of lowest MSFE'};
T13 = table(SBL,SNS,LASSO,SSVS,BMA,BAG,DFM5,OLS,'RowNames',Name);