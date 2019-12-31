% Forecasting using Bayesian Hierarchical Shrinkage regressions
%==========================================================================
% This code replicates the forecasting results in Table 3 of Korobilis (2019)
% "High-dimensional macroeconomic forecasting using message passing
% algorithms"
% 
% Instructions: 
% - Make sure you choose CPI when running this code (set "index_variable = 95" in USER INPUT).
% - You need to run the code for each forecast horizon (set "nfore = h" where h is the desired horizon 1,3,6,12)
% - If code runs too slow, try commenting the TVP_TVD and TVP_TVS functions (in Estimation section below, see 1d 
%   around line 145) and all subsequent lines of code that use output from these functions (e.g. y_hat_TVPTVD and 
%   y_hat_TVPTVS etc). Everything else should run reasonably fast.
%==========================================================================
 
% Reset everything, clear memory and screen
clear; close all; clc;
% Start clock
tic;
% Add path of random number generators
addpath('functions')
addpath('data')
 
%===========================| USER INPUT |=================================
% Model specification
constant = 1;           % 0: no intercept; 1: intercept
lags = 2;               % Number of AR lags
 
% Forecasting
forecasting = 1;        % 0: no forecasting; 1: forecasting
index_variable = 95;    % 95:CPI; 105:PCE;
nfore = 12;             % Number of predictions steps ahead
fore_type = 0;          % 0: Recursive, 1: Rolling
ndraws = 5000;          % Number of draws
nburn  = 1000;          % Burn in period for MCMC draws
%% ===========================| LOAD DATA |=================================
% Read in the data
% [A,B,~] = xlsread('FREDQD.xlsx','FREDQD');
% data    = A(3:end,:);
% tcode   = A(2,:);      % Vector of transformations
% fcode   = A(1,:);      % Variables to extract PCA from
% Ynames  = B(1,2:end);  % The mnemonic for the variables in the panel
% imp_pred = [1,219,140,120,161,194,32,73,174,72] + constant + lags; %10 important predictors
 
[A,B,~] = xlsread('FREDMD.xlsx','Data');
data    = A(2:end,:);
tcode   = A(1,:);      % Vector of transformations
Ynames  = B(1,2:end);  % The mnemonic for the variables in the panel
imp_pred = [6,64,82,93,87,114,19,48,53,20] + constant + lags; %10 important predictors
%---------------
 
% Start with  first_per% of observations for first period forecast
first_per = 0.5;
T_full    = (size(data,1)-2-nfore) - nfore;
T_thres   = round(first_per*T_full); 
anumber   = T_full-T_thres+1;
 
MSFE  = zeros(anumber,10);
MAFE  = zeros(anumber,10);
logPL = zeros(anumber,10);
 
% Choose one of the series above
Yraw_full = yfcsta(transxf(data(:,index_variable),5),5,nfore);   Yraw_full = adjout(Yraw_full,4.5,4);
Yraw_full_lags = transxf(data(:,index_variable),5);              Yraw_full_lags = adjout(Yraw_full_lags,4.5,4);
Xraw = 0*data;
for kk = 1:size(data,2)
    Xraw(:,kk) = transxf(data(:,kk),tcode(kk));
    Xraw(:,kk) = adjout(Xraw(:,kk),4.5,4);
end
% a) Remove the dependent variables from the set of predictors (and correct location of 10 important predictors of inflation)
Xraw(:,index_variable) = []; imp_pred((imp_pred-constant-lags)>index_variable)=imp_pred((imp_pred-constant-lags)>index_variable)-1;
% b) Correct for transformation (for second log-diffs we lose 2 observations)
Xraw = Xraw(3:end,:);  Yraw_full = Yraw_full(3:end,:); Yraw_full_lags = Yraw_full_lags(3:end,:);
 
% Now correct sizes for h-step ahead forecasts
Yraw_full = Yraw_full(1:end-nfore,:);  
Yraw_full_lags = Yraw_full_lags(1:end-nfore,:);
Xraw = Xraw(1:end-nfore,:);
 
for sample = T_thres:T_full
    disp('     ')
    disp(['Now you are running sample ' num2str(sample) ' of ' num2str(T_full)] )
    disp(['Variable you are using is ' Ynames(index_variable) 'which is number ' num2str(index_variable)])           
     
    if fore_type == 0 % Recursive forecasts
        roll = 1;
    elseif fore_type == 1 % Rolling forecasts         
        roll = sample-99; 
    end
     
    p = 20;%size(Xraw,2);   
    q = constant+lags;
    T = sample;
    % Create RHS variables, x
    ylag = mlag2(Yraw_full_lags,lags-1);  % This routine creates the lagged dep variables
%     F = extract(zscore(Xraw(1:T+nfore,find(fcode==1))),p);
    F = extract(zscore(Xraw(1:T+nfore,:)),p);
    %F = pc_factor(zscore(Xraw(1:T+nfore,:)),p);
    if constant == 1       
        Zraw  = [ones(size(Yraw_full_lags(lags:sample+nfore,:),1),1) Yraw_full_lags(lags:sample+nfore,:) ylag(lags:sample+nfore,:) Xraw(lags:sample+nfore,:)];
        Zraw2 = [ones(size(Yraw_full_lags(lags:sample+nfore,:),1),1) Yraw_full_lags(lags:sample+nfore,:) ylag(lags:sample+nfore,:) F(lags:sample+nfore,:)];
    elseif constant == 0                        
        Zraw  = [Yraw_full_lags(lags:sample+nfore,:) ylag(lags:sample+nfore,:) Xraw(lags:sample+nfore,:)];
        Zraw2 = [Yraw_full_lags(lags:sample+nfore,:) ylag(lags:sample+nfore,:) F(lags:sample+nfore,:)];
    end
    % Correct # of observations of y after taking lags
    Yraw = Yraw_full(lags:sample,:);
     
    %In-sample observations
    y     = Yraw(roll:end,:); ss = std(y); mm = mean(y); y = zscore(y);
    Zraw(:,1+constant:end) = zscore(Zraw(:,1+constant:end)); Zraw2(:,1+constant:end) = zscore(Zraw2(:,1+constant:end));
%     columnNorms = sqrt(diag(Zraw'*Zraw)); Zraw = Zraw*diag(1 ./ columnNorms);
    Xtemp = Zraw(roll:end,[1:q, imp_pred]);  Xtemp2 = Zraw2(roll:end,:);
    X     = Xtemp(1:end-nfore,:); X2    = Xtemp2(1:end-nfore,:);
    T     = size(y,1);
       
    %Out-of-sample observations
    Yraw_fore = Yraw_full(sample+nfore,:);
    X_fore    = Xtemp(end,:);  X_fore2   = Xtemp2(end,:);
     
    for i = 1:T
        XX(i,(i-1)*(p+q)+1:i*(p+q))    = X2(i,:);
        XX2(i,(i-1)*q+1:i*q)           = X(i,1:q);
        XX3(i,(i-1)*(q+10)+1:i*(q+10)) = X(i,:);
        XX_fore(1,(i-1)*(p+q)+1:i*(p+q))    = X_fore2;
        XX_fore2(1,(i-1)*q+1:i*q)           = X_fore(:,1:q);   
        XX_fore3(1,(i-1)*(q+10)+1:i*(q+10)) = X_fore;
    end
     
    %% Estimation
    % 1) Models with no predictors (just AR terms)
    % a\ constant parameter AR
    [beta_ARp,sigma_ARp]             = OLS(X(:,1:q),y);
    
    % b\ structural breaks AR
    [y_f_KP]                         = KP_PPT(y,X(:,1:q),X_fore(:,1:q),1,ndraws,nburn);
    [y_f_GK]                         = TVP_GK(y,X(:,1:q),X_fore(:,1:q),1,ndraws,nburn);
    
    % c\ tvp AR
    [beta_TVPAR,sigma_TVPAR]         = TVP_AR(y,X(:,1:q),ndraws,nburn);
    [beta_UCSV,sigma_UCSV]           = TVP_UCSV(y,ndraws,nburn);       
    
    % d\ tvp AR with shrinkage
    [y_hat_TVPTVD]                   = TVP_TVD(y,X2(:,1:q+3),X_fore2(:,1:q+3),0,ndraws/2,nburn);  % Chan et al (2012) JBES
    [beta_TVPTVS,sigma_TVPTVS]       = TVP_TVS(y,X2(:,1:q+3),0.1,0.1,ndraws/2,nburn);        % Kalli+Griffin (2014) JoE    
    
    % 2) Models with predictors
    % a\ constant parameter regressions
    [beta_SSVS,~,sigma_SSVS,~]       = SSVS(X2,y,ndraws,nburn);   % This is called BMA in the paper (because SSVS can be used either for BMA or BMS)

    % b\ tvp regressions
    [beta_TVPBMA,sigma_TVPBMA]       = TVP_BMA(y,X2(:,1:13),lags,ndraws,nburn);
    
    % c\ TVP-GAMP
    [beta_SBLTV,~,sigma_SBLTV,~]     = AMPSBL_SV([X2, XX],y,1000,1e-10,1e-10);
     
                   
    %% Forecasting
    if forecasting == 1
        y_hat_ARp     = mm + ss.*(X_fore(:,1:q)*beta_ARp + sqrt(sigma_ARp)*randn(1,ndraws));
        y_hat_KP      = mm + ss.*y_f_KP';
        y_hat_GK      = mm + ss.*y_f_GK';
        y_hat_TVPAR   = mm + ss.*(X_fore(:,1:q)*squeeze(beta_TVPAR(end,:,:)) + sqrt(sigma_TVPAR(end,:)).*randn(1,ndraws));
        y_hat_UCSV    = mm + ss.*(beta_UCSV(end,:) + sqrt(sigma_UCSV(end,:)).*randn(1,ndraws));
        y_hat_TVPTVD  = mm + ss.*(y_hat_TVPTVD');
        y_hat_TVPTVS  = mm + ss.*(X_fore2(:,1:q+3)*squeeze(beta_TVPTVS(end,:,:)) + sqrt(sigma_TVPTVS(end,:)).*randn(1,ndraws/2));
        y_hat_SSVS    = mm + ss.*(X_fore2*beta_SSVS + sqrt(sigma_SSVS)*randn(1,ndraws));
        y_hat_TVPBMA  = mm + ss.*(X_fore2(:,1:13)*squeeze(beta_TVPBMA(end,:,:)) + sqrt(sigma_TVPBMA(end,:)).*randn(1,ndraws));
        y_hat_SBLTV   = mm + ss.*([X_fore2, XX_fore]*beta_SBLTV + sqrt(sigma_SBLTV(end))*randn(1,ndraws));
         
        PL_ARp        = ksdensity(y_hat_ARp,Yraw_fore);
        PL_KP         = ksdensity(y_hat_KP,Yraw_fore);
        PL_GK         = ksdensity(y_hat_GK,Yraw_fore);
        PL_TVPAR      = ksdensity(y_hat_TVPAR,Yraw_fore);
        PL_UCSV       = ksdensity(y_hat_UCSV,Yraw_fore);
        PL_TVPTVD     = ksdensity(y_hat_TVPTVD,Yraw_fore);
        PL_TVPTVS     = ksdensity(y_hat_TVPTVS,Yraw_fore);
        PL_SSVS       = ksdensity(y_hat_SSVS,Yraw_fore);
        PL_TVPBMA     = ksdensity(y_hat_TVPBMA,Yraw_fore);
        PL_SBLTV      = ksdensity(y_hat_SBLTV,Yraw_fore);
                
        med_f_ARp     = mean(y_hat_ARp',1);       
        med_f_KP      = mean(y_hat_KP',1);
        med_f_GK      = mean(y_hat_GK',1);
        med_f_TVPAR   = mean(y_hat_TVPAR',1);
        med_f_UCSV    = mean(y_hat_UCSV',1);
        med_f_TVPTVD  = median(y_hat_TVPTVD',1);
        med_f_TVPTVS  = median(y_hat_TVPTVS',1);
        med_f_SSVS    = mean(y_hat_SSVS',1);
        med_f_TVPBMA  = mean(y_hat_TVPBMA',1);
        med_f_SBLTV   = median(y_hat_SBLTV',1);
         
        % Now calculate the statistics of interest
        MSFE(sample-T_thres+1,1)   = (med_f_ARp - Yraw_fore).^2;
        MSFE(sample-T_thres+1,2)   = (med_f_KP - Yraw_fore).^2;
        MSFE(sample-T_thres+1,3)   = (med_f_GK - Yraw_fore).^2;
        MSFE(sample-T_thres+1,4)   = (med_f_TVPAR - Yraw_fore).^2;      
        MSFE(sample-T_thres+1,5)   = (med_f_UCSV - Yraw_fore).^2;
        MSFE(sample-T_thres+1,6)   = (med_f_TVPTVD - Yraw_fore).^2;
        MSFE(sample-T_thres+1,7)   = (med_f_TVPTVS - Yraw_fore).^2;
        MSFE(sample-T_thres+1,8)   = (med_f_SSVS - Yraw_fore).^2;
        MSFE(sample-T_thres+1,9)   = (med_f_TVPBMA - Yraw_fore).^2;
        MSFE(sample-T_thres+1,10)  = (med_f_SBLTV - Yraw_fore).^2;
         
        MAFE(sample-T_thres+1,1)   = abs(med_f_ARp - Yraw_fore);       
        MAFE(sample-T_thres+1,2)   = abs(med_f_KP - Yraw_fore);
        MAFE(sample-T_thres+1,3)   = abs(med_f_GK - Yraw_fore);
        MAFE(sample-T_thres+1,4)   = abs(med_f_TVPAR - Yraw_fore);
        MAFE(sample-T_thres+1,5)   = abs(med_f_UCSV - Yraw_fore);
        MAFE(sample-T_thres+1,6)   = abs(med_f_TVPTVD - Yraw_fore);
        MAFE(sample-T_thres+1,7)   = abs(med_f_TVPTVS - Yraw_fore);
        MAFE(sample-T_thres+1,8)   = abs(med_f_SSVS - Yraw_fore);
        MAFE(sample-T_thres+1,9)   = abs(med_f_TVPBMA - Yraw_fore);
        MAFE(sample-T_thres+1,10)  = abs(med_f_SBLTV - Yraw_fore);        
         
        logPL(sample-T_thres+1,1)  = mean(log(PL_ARp + 1e-16));
        logPL(sample-T_thres+1,2)  = mean(log(PL_KP + 1e-16));
        logPL(sample-T_thres+1,3)  = mean(log(PL_GK + 1e-16));
        logPL(sample-T_thres+1,4)  = mean(log(PL_TVPAR + 1e-16));
        logPL(sample-T_thres+1,5)  = mean(log(PL_UCSV + 1e-16));
        logPL(sample-T_thres+1,6)  = mean(log(PL_TVPTVD + 1e-16));
        logPL(sample-T_thres+1,7)  = mean(log(PL_TVPTVS + 1e-16));
        logPL(sample-T_thres+1,8)  = mean(log(PL_SSVS + 1e-16));
        logPL(sample-T_thres+1,9)  = mean(log(PL_TVPBMA + 1e-16));
        logPL(sample-T_thres+1,10) = mean(log(PL_SBLTV + 1e-16));
    end
end %recursive forecasts

save(sprintf('%s_%s_%g_%s_%g.mat',cell2mat(Ynames(index_variable)),'(_h_=',nfore,')',rand),'-mat');