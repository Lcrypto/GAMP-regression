function [beta_means, sigma2] = BAGGING(x,y,n_boots)

Tcrit = 2.807;
bmat= [x,y];
npred=size(bmat,2)-1; %number of predictors           
T2=size(bmat,1);      %number of observations %here
SelectP=cell(n_boots,1);
Predmat=cell(n_boots,1);           
forec_b=zeros(n_boots,1);
resid_boot=cell(1,n_boots);         
var_boot=cell(n_boots,1);
betas_all = zeros(npred,n_boots);
sigmas_all = zeros(n_boots,1);
for j=1:n_boots
    if j==1  
        DS=bmat;
    else
        %sampling should be done in blocs
        % ceil(T/m) Uniform random numbers over 1...T-m+1
        Index = ceil((T2)*rand(ceil(T2),1));
        Index = bsxfun(@plus,Index,0)';
        % Transform to col vector, and remove excess
        Index = Index(:); Index = Index(1:T2);
        DS= bmat(Index,(1:npred+1));  %Construct the data
    end    
    Predmat{j,1}=DS;
           
    %Estimate the multiple regression with HAC standard errors
    [Coef_bt,StE] = myols(DS(1:T2,npred+1),[ones(T2,1),DS(1:T2,1:npred)]);
    %Compute the T-stat
    Tstat=Coef_bt./StE;
    
    %select variables with T>|Tcrit|
    SelectP{j,1}= find((abs(Tstat(2:npred+1)))>Tcrit);
    
    %Forecasting by picking the predictors selected in mimmodel
    varsel=SelectP{j,1};
    DS1=Predmat{j,1};
                
    if isempty(varsel)==1
        forec_b(j,1)=0;
        resid_boot{1,j}=DS1(:,npred+1);
        Coef_fc=0;
        sigma2_j = (DS1(1:T2,npred+1))'*(DS1(1:T2,npred+1))/(T2-numel(varsel)-1);
    else        
        Newbmat=[ones(T2,1), DS1(1:T2,varsel)];
        Coef_fc=Newbmat\DS1(1:T2,npred+1);
        sigma2_j = (DS1(1:T2,npred+1) - Newbmat*Coef_fc)'*(DS1(1:T2,npred+1) - Newbmat*Coef_fc)/(T2-numel(varsel)-1);
    end
    betas_all([1 ; varsel+1],j) = Coef_fc;
    sigmas_all(j,1) = sigma2_j;
end

beta_means = mean(betas_all,2);
sigma2 = mean(sigmas_all);
