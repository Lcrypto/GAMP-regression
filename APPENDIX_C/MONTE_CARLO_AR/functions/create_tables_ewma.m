function [T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13] = create_tables_ewma(sel_series)

if nargin == 0
    sel_series = 1:222;
end


% MSFEs
% Table 1
load('C:\Users\Dell\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\ewma_ma\FORECASTING_EWMA_(_h_=_1_).mat')
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,4,sel_series))),1,4));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MRMSFE([1:3],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); DFM5 = entries(:,3);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T1 = table(SBL,SNS,DFM5,'RowNames',Name);

% Table 2
load('C:\Users\Dell\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\ewma_ma\FORECASTING_EWMA_(_h_=_2_).mat')
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,4,sel_series))),1,4));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MRMSFE([1:3],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); DFM5 = entries(:,3);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T2 = table(SBL,SNS,DFM5,'RowNames',Name);

% Table 3
load('C:\Users\Dell\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\ewma_ma\FORECASTING_EWMA_(_h_=_4_).mat')
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,4,sel_series))),1,4));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MRMSFE([1:3],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); DFM5 = entries(:,3);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T3 = table(SBL,SNS,DFM5,'RowNames',Name);

% Table 4
load('C:\Users\Dell\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\ewma_ma\FORECASTING_EWMA_(_h_=_8_).mat')
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,4,sel_series))),1,4));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MRMSFE([1:3],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); DFM5 = entries(:,3);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T4 = table(SBL,SNS,DFM5,'RowNames',Name);

% MAFEs
% Table 5
load('C:\Users\Dell\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\ewma_ma\FORECASTING_EWMA_(_h_=_1_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,4,sel_series)),1,4));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MMAFE([1:3],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); DFM5 = entries(:,3);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T5 = table(SBL,SNS,DFM5,'RowNames',Name);

% Table 6
load('C:\Users\Dell\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\ewma_ma\FORECASTING_EWMA_(_h_=_2_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,4,sel_series)),1,4));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MMAFE([1:3],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); DFM5 = entries(:,3);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T6 = table(SBL,SNS,DFM5,'RowNames',Name);

% Table 7
load('C:\Users\Dell\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\ewma_ma\FORECASTING_EWMA_(_h_=_4_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,4,sel_series)),1,4));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MMAFE([1:3],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); DFM5 = entries(:,3);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T7 = table(SBL,SNS,DFM5,'RowNames',Name);

% Table 8
load('C:\Users\Dell\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\ewma_ma\FORECASTING_EWMA_(_h_=_8_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,4,sel_series)),1,4));  % Mean Relative MSFE (benchmark is ARp)
entries = quantile(MMAFE([1:3],:)',[0.05 0.25 0.50 0.75 0.95]);
SBL = entries(:,1); SNS = entries(:,2); DFM5 = entries(:,3);
Name = {'0.05', '0.25', '0.50', '0.75', '0.95'};
T8 = table(SBL,SNS,DFM5,'RowNames',Name);

% Correlations
% Table 9
load('C:\Users\Dell\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\ewma_ma\FORECASTING_EWMA_(_h_=_1_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,4,sel_series)),1,4));
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,4,sel_series))),1,4));  % Mean Relative MSFE (benchmark is ARp)
RMSFE = zeros(3,3);
for i=1:2
    for j=i+1:3
        RMSFE(i,j) = mean(abs(MRMSFE(i,:)-MRMSFE(j,:)));
    end
end
entries = tril(corrcoef(MMAFE(1:3,:)')) + RMSFE - eye(3);
SBL = entries(:,1); SNS = entries(:,2); DFM5 = entries(:,3);
Name = {'SBL','SNS','DFM5'};
T9 = table(SBL,SNS,DFM5,'RowNames',Name);

% Hit Rates
% Table 10
load('C:\Users\Dell\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\ewma_ma\FORECASTING_EWMA_(_h_=_1_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,4,sel_series)),1,4));  % Mean Relative MSFE (benchmark is ARp)
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,4,sel_series))),1,4));  % Mean Relative MSFE (benchmark is ARp)
[~,b1] = min(MMAFE(1:3,:)); [~,b2] = min(MRMSFE(1:3,:));
SBL = 100*[sum(b1==1)/length(sel_series);sum(b2==1)/length(sel_series)] ; SNS = 100*[sum(b1==2)/length(sel_series);sum(b2==2)/length(sel_series)]; DFM5 = 100*[sum(b1==3)/length(sel_series);sum(b2==3)/length(sel_series)];
Name = {'% of lowest MAFE','% of lowest MSFE'};
T10 = table(SBL,SNS,DFM5,'RowNames',Name);

% Table 11
load('C:\Users\Dell\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\ewma_ma\FORECASTING_EWMA_(_h_=_2_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,4,sel_series)),1,4));  % Mean Relative MSFE (benchmark is ARp)
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,4,sel_series))),1,4));  % Mean Relative MSFE (benchmark is ARp)
[~,b1] = min(MMAFE(1:3,:)); [~,b2] = min(MRMSFE(1:3,:));
SBL = 100*[sum(b1==1)/length(sel_series);sum(b2==1)/length(sel_series)] ; SNS = 100*[sum(b1==2)/length(sel_series);sum(b2==2)/length(sel_series)]; DFM5 = 100*[sum(b1==3)/length(sel_series);sum(b2==3)/length(sel_series)];
Name = {'% of lowest MAFE','% of lowest MSFE'};
T11 = table(SBL,SNS,DFM5,'RowNames',Name);

% Table 12
load('C:\Users\Dell\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\ewma_ma\FORECASTING_EWMA_(_h_=_4_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,4,sel_series)),1,4));  % Mean Relative MSFE (benchmark is ARp)
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,4,sel_series))),1,4));  % Mean Relative MSFE (benchmark is ARp)
[~,b1] = min(MMAFE(1:3,:)); [~,b2] = min(MRMSFE(1:3,:));
SBL = 100*[sum(b1==1)/length(sel_series);sum(b2==1)/length(sel_series)] ; SNS = 100*[sum(b1==2)/length(sel_series);sum(b2==2)/length(sel_series)]; DFM5 = 100*[sum(b1==3)/length(sel_series);sum(b2==3)/length(sel_series)];
Name = {'% of lowest MAFE','% of lowest MSFE'};
T12 = table(SBL,SNS,DFM5,'RowNames',Name);

% Table 13
load('C:\Users\Dell\Dropbox\Research\GAMP\Code\2017.04.30\OUTPUT\FORECASTING\ewma_ma\FORECASTING_EWMA_(_h_=_8_).mat')
MMAFE = squeeze(mean(MAFE(:,:,sel_series))./repmat(mean(MAFE(:,4,sel_series)),1,4));  % Mean Relative MSFE (benchmark is ARp)
MRMSFE = squeeze(sqrt(mean(MSFE(:,:,sel_series)))./repmat(sqrt(mean(MSFE(:,4,sel_series))),1,4));  % Mean Relative MSFE (benchmark is ARp)
[~,b1] = min(MMAFE(1:3,:)); [~,b2] = min(MRMSFE(1:3,:));
SBL = 100*[sum(b1==1)/length(sel_series);sum(b2==1)/length(sel_series)] ; SNS = 100*[sum(b1==2)/length(sel_series);sum(b2==2)/length(sel_series)]; DFM5 = 100*[sum(b1==3)/length(sel_series);sum(b2==3)/length(sel_series)];
Name = {'% of lowest MAFE','% of lowest MSFE'};
T13 = table(SBL,SNS,DFM5,'RowNames',Name);