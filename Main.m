%HAYASHI EXERCISE 1: b,d,e,f,g,h,i,j,k,l
%Mishkin has been included in this folder
%JohnDaniel Paletto
%Last modified 10/21/15






%!@#$%^&*()PLEASE REFRAIN FROM CLICKING "RUN", PLEASE USE "RUN SECTION"!@#$%^&*()
%%%%% THIS WILL RUN PART BY PART and avoid DEATH-BY-POP-UPS






clear       %delete/clear memory
clc         %clear output screen
close all   %close e.g. figures
load('mishkin.mat');
%Notes from Mishkin.mat: By Column: Year,Month, 1mo inflation rate, 3mo
%inflation rate, 1mo Tbill, 3mo Tbill, CPI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PART (b)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The CPI is consistent with Urban Consumer data (all items) from BLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END PART (b)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%START PART (d)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First we will calculate a 1mo inflation rate from Mishkin CPI
for i=1:490                 
    cpir(i+1,1)=((Mishkin(i+1,7)/Mishkin(i,7))^12-1)*100;
end;

%Calculate real rate by pairing 1mo inflation rate above to 1mo Tbill
tb1=Mishkin(:,5);
r=tb1-cpir;

%We want to Duplicate FAMA(1975), So truncate & Assume we know (t-1)
r=r(36:258,1);

%Calculate auto Correlation using autocov and autocorell(x,maxlag)
pj=autocorrel(r,12);

%Calculate STD ER of the Pj's
stde=((1-pj.^2)/(length(r))).^(1/2);

%LJUNG Q: (vector x, maxlag), returns vector of q statistics per lag
q=LjungQ(r, 12);

%CHi PVAL for each lag
for l=1:12
p(l)=(1-chi2cdf(q(l),l))*100;
end;

%PLOTS
x=1:491;
ze=zeros(491);

tbilvsCPIrate=figure;               %tbill 1mo rate red
set(tbilvsCPIrate, 'name', '1mo Tbill rate-RED / CPI inflation rate-BLUE')
tbilvsCPIrate,
plot(x,cpir,x,tb1,x,ze);            %CPI inflation rate blue
hold off

x2=1:223;
ze2=zeros(223);

Realrate=figure;
set(Realrate, 'name', 'Ex-post Real Rate')
Realrate,
plot(x2,r,x2,ze2);
hold off

%Figure replication of Table 2.1
rirft(1,:)=pj';
rirft(2,:)=stde';
rirft(3,:)=q';
rirft(4,:)=p;

rnames = {'phatj','Std Error','LJung Q','P-Val'}; 
cnames = {'j=1','j=2','j=3','j=4','j=5','j=6','j=7','j=8','j=9','j=10',...
    'j=11','j=12'}

f=figure('Position', [100 100 1000 200]);
set(f, 'name', 'Real Interest Rates, January 1953-July 1971',...
    'numbertitle','off');

realInterestRateFamaTable = uitable('Data',rirft,'RowName',rnames,...
    'ColumnName',cnames,...
    'Tag','Real Interest Rates, January 1953-July 1971',...
    'TooltipString','Real Interest Rates, January 1953-July 1971',...
    'Parent', f,'Position',[50 50 1000 150]);

realInterestRateFamaTable

%NOTES
%These tables/process/vaalues repoduce the exact values found in the text
%We can therefore affirm our process is 'most likely' correct

%%%%%%%%%%%%%%%%%%%%%%%%%END of PART (d)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%START PART (e)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We want to run a regression where x = (1,tb1)' and y = cpir
X(:,1)=ones(223,1);      %intercept
X(:,2)=tb1(36:258);      %t-1 tbill 1mo rate
y=cpir(36:258);          %t CPI inflation rate

beta_hat = ((X'*X)^(-1))*(X'*y);                     %regression coefs
e_hat = y-X*beta_hat;                                %estimation of resids
mean_y=mean(y);                                      %calc mean y
R_sqr = 1-(e_hat'*e_hat)/((y-mean_y)'*(y-mean_y));   %R sqr for regression
n=length(y);
B = diag(e_hat.^2);                                  %Create diag of res^2
S_hat = (X'*B*X)/(n-1);                              %sum xixi'*e^2/n
Sxx = X'*X/(n-1);                                    %Variance matrix of X
Avar_hat_betas = Sxx^(-1)*S_hat*Sxx^(-1);            %Asymtotic variance X

SE_B0 = (Avar_hat_betas(1,1)/(n-1))^(1/2);               %Std error of intcpt
SE_B1 = (Avar_hat_betas(2,2)/(n-1))^(1/2);               %Std error of Beta tb1
SE_R = ((e_hat-mean(e_hat))'*(e_hat-mean(e_hat))/(n-2))^(1/2) %Std error Re

%ROBUST T-STAT
robust_t_beta1_OLS = beta_hat(2,1)/SE_B1

%Optional T-test
%   Process By: Daniela Osterrieder, Rutgers University, 2015; 
%   Only minor changes were made and commnents left intact
%H0: beta = 0 vs H1: beta>0 
critical_value = icdf('normal',0.95,0,1); 
if robust_t_beta1_OLS<critical_value
    'based on a significance test with robust SEs, we fail to reject beta=0'
    
else
    'based on a significance test with robust SEs, we reject beta=0 in favour of the alternative hypothesis that beta>0'
end;
%   end Osterrieder, 2015

%PLOTS & Tables
reg=figure;
set(reg,'name','Scatter and Best fit line FAMA REGRESSION TB1')
reg,
scatter(X(:,2),y);
hold on
plot(X(:,2),1.0147*X(:,2)-.8678);

R_tbl(1,1:2)=beta_hat;
R_tbl(2,1:2)=[SE_B0,SE_B1];
R_tbl(3,3:6)=[R_sqr,mean_y,SE_R,robust_t_beta1_OLS];

freg=figure('Position', [150 150 600 255]);
set(freg, 'name', 'FAMA Regression Stats w/ TB1','numbertitle','off');

r2names={'Coefficients','Std Error','Statistics'};
c2names={,'B0','B1','R^2','Mean of y', 'SER in %','T-stat'};
Regression_table = uitable('Data',R_tbl,...
    'RowName',r2names,'ColumnName',c2names,'Tag',...
    'Fama Regression Inflation Rt by Tbill Rt, January 1953-July 1971',...
    'Parent', freg,'Position',[40 40 550 155]);

%NOTES
%The intercept accounts for all excluded variables, and for that of which
%is 'unexplained' in the regression.
%
%The table and regression reproduced are exactly the same as the text. We
%have used an estimate for the standard errors of rho. These values are
%very close to the text.

%%%%%%%%%%%%%%%%%%%%%%%%%%END OF PART(e)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%START of PART (f)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We want to further evaluate the standard errors
%FIRST, we can use a degrees of freedom correction *n/(n-K)
df_corrected_S_hat = S_hat*(n-1)/(n-2);
temp = Sxx^(-1)*df_corrected_S_hat*Sxx^(-1);
df_corrected_SE_B0 = (temp(1,1)/(n-1))^(1/2);
df_corrected_SE_B1 = (temp(2,2)/(n-1))^(1/2);
clear temp

%SECOND, we can use 2.5.5 with d = 1
for i = 1:223                                        %Loop to calc p diag
    piv01(i,1) = X(i,1:2)*(X'*X)^(-1)*X(i,1:2)';
end;

d1_DMac_S_hat=zeros(2,2);
for i = 1:223
    d1_DMac_S_hat=d1_DMac_S_hat+((e_hat(i)^2/(1-piv01(i)))*X(i,1:2)'*X(i,1:2));
end;
d1_DMac_S_hat=d1_DMac_S_hat./n;
temp = Sxx^(-1)*d1_DMac_S_hat*Sxx^(-1);

d1_corrected_SE_B0 = (temp(1,1)/n)^(1/2);
d1_corrected_SE_B1 = (temp(2,2)/n)^(1/2);
clear temp

%THIRD, we can use 2.5.5 with d=2
d2_DMac_S_hat=zeros(2,2);
for i = 1:223
    d2_DMac_S_hat=d2_DMac_S_hat+((e_hat(i)^2/(1-piv01(i)^2))*X(i,1:2)'*X(i,1:2));
end
d2_DMac_S_hat=d2_DMac_S_hat./n;
temp = Sxx^(-1)*d2_DMac_S_hat*Sxx^(-1);

d2_corrected_SE_B0 = (temp(1,1)/n)^(1/2);
d2_corrected_SE_B1 = (temp(2,2)/n)^(1/2);
clear temp

%NOTES
%These corrected SE are very close to the actual which we would assume to be true.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END PART (f)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%START PART (g)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Regression under the assumption of conditional Homoskedasticity
S_hat_HOMO = SE_R^2*Sxx;
Avar_b_hat_HOMO = n*SE_R^2*(X'*X)^(-1);              %Asymtotic variance X
SE_B0_HOMO = (Avar_b_hat_HOMO(1,1)/n)^(1/2);         %Std error of intcpt
SE_B1_HOMO = (Avar_b_hat_HOMO(2,2)/n)^(1/2)         %Std error of Beta tb1

%T-RATIO FROM 1.4.5
t_beta1_OLS_HOMO = beta_hat(2,1)/SE_B1_HOMO

%Optional T-test
%   Process By: Daniela Osterrieder, Rutgers University, 2015; 
%   Only minor changes were made and commnents left intact
%H0: beta = 0 vs H1: beta>0 
critical_value = icdf('normal',0.95,0,1); 
if t_beta1_OLS_HOMO<critical_value
    'based on a significance test with non-robust SEs, we fail to reject beta=0'
    
else
    'based on a significance test with non-robust SEs, we reject beta=0 in favour of the alternative hypothesis that beta>0'
end;
%   end Osterrieder, 2015

R_tblhomo(1,1:2)=beta_hat;
R_tblhomo(2,1:2)=[SE_B0_HOMO,SE_B1_HOMO];
R_tblhomo(3,3:6)=[R_sqr,mean_y,SE_R,t_beta1_OLS_HOMO];

freghomo=figure('Position', [150 150 600 255]);
set(freghomo, 'name', 'FAMA Regression Stats w/ TB1 under Homoskedasticity','numbertitle','off');

r2names={'Coefficients','Std Error','Statistics'};
c2names={,'B0','B1','R^2','Mean of y', 'SER in %','T-stat'};
Regression_table = uitable('Data',R_tblhomo,...
    'RowName',r2names,'ColumnName',c2names,'Tag',...
    'Fama Regression Inflation Rt by Tbill Rt, January 1953-July 1971',...
    'Parent', freghomo,'Position',[40 40 550 155]);

%NOTES
%Estimators are the same.
%S-hat is different therefor variance is different along with the standard
%errors and t-statistic. Although the regressions are similar, the robust 
%t is larger

%%%%%%%%%%%%%%%%%%%%%%%%%%%END PART (h)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We will conduct a Breusch-Godfrey

%   Process By: Daniela Osterrieder, Rutgers University, 2015; ends ~35 lines
%   Only minor changes were made and commnents left intact
bg=[1:1:223]';
Breusch_Godfrey=NaN(length(pj),1);
for j=12:-1:1
    maxlag=bg(j,1);
    independent=X;
    for i=1:maxlag
        temp=[zeros(i,1); e_hat(1:end-i,1)]; %set presample values equal to zero
        independent=[independent temp];
        clear temp
    end;
    
    dependent=e_hat;
    
    coeff_auxil = ((independent'*independent)^(-1))*(independent'*dependent);
    resid_auxil = dependent-independent*coeff_auxil;
    clear coeff_auxil independent
    
    R2=1-((resid_auxil'*resid_auxil)/((dependent-mean(dependent))'*(dependent-mean(dependent))));
    clear dependent
    Breusch_Godfrey(j,1)=(n)*R2;
    
    critical_value_Q(j,1) = icdf('chi2',0.95,maxlag);
    if Breusch_Godfrey(j,1)<critical_value_Q(j,1)
        'based on the Breusch-Godfrey test for serial correlation up to lag' 
        maxlag 
        ', we fail to reject the assumption of no serial correlation for regression from part (d)'
        '-----------------------------------------------------------------'
    else
        'based on the Breusch-Godfrey test for serial correlation up to lag' 
        maxlag 
        ', we reject the assumption of no serial correlation for regression from part (d)'
        '-----------------------------------------------------------------'
        break
        %If you desire all 12 BG statistics remove the break
    end;
end;
%     end of Osterrieder, 2015

%PLOTS & TABLES
figure,
bar(pj)
hold off

%NOTES
%BG Statistic
R2*n
%27, which is the value mentioned in Hayashi
%So, we can affirm that our process and code produce the correct statistic
%at least for p=12

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END PART (h)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%START PART (j)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We add monthly seasonal dummies instead of a constant
dummies=zeros(223,12);                                %initialize
dummies=repmat(eye(12,12),18,1);                      %create dummy matrix
dummies(217:223,1:7)=eye(7,7);

X_w_dummies(1:223,1)=X(:,2);                          %create X matrix
X_w_dummies(1:223,2:13)=dummies;                      %add dummies

dum_beta_hat = ((X_w_dummies'*X_w_dummies)^(-1))*(X_w_dummies'*y);

e_hat_dum = y-X_w_dummies*dum_beta_hat;                   %estimation of resids
R_sqr_dum = 1-(e_hat_dum'*e_hat_dum)/((y-mean_y)'*(y-mean_y));   %fit using R^2

B_dum = diag(e_hat_dum.^2);                           %Create diag of res^2
S_hat_dum = (X_w_dummies'*B_dum*X_w_dummies)/(n-13);       %sum xixi'*e^2/n
Sxx_dum = X_w_dummies'*X_w_dummies/(n-13);            %Variance matrix of X
Avar_hat_b_dum = Sxx_dum^(-1)*S_hat_dum*Sxx_dum^(-1);   %Asymtotic variance X

SE_tb1_dum = (Avar_hat_b_dum(1,1)/(n-1))^(1/2);        %Std error of tb1
SE_M1_dum = (Avar_hat_b_dum(2,2)/(n-1))^(1/2);         %Std error of M1
SE_M2_dum = (Avar_hat_b_dum(3,3)/(n-1))^(1/2);         %Std error of M2
SE_M3_dum = (Avar_hat_b_dum(4,4)/(n-1))^(1/2);         %Std error of M3
SE_M4_dum = (Avar_hat_b_dum(5,5)/(n-1))^(1/2);         %Std error of M4
SE_M5_dum = (Avar_hat_b_dum(6,6)/(n-1))^(1/2);         %Std error of M5
SE_M6_dum = (Avar_hat_b_dum(7,7)/(n-1))^(1/2);         %Std error of M6
SE_M7_dum = (Avar_hat_b_dum(8,8)/(n-1))^(1/2);         %Std error of M7
SE_M8_dum = (Avar_hat_b_dum(9,9)/(n-1))^(1/2);         %Std error of M8
SE_M9_dum = (Avar_hat_b_dum(10,10)/(n-1))^(1/2);       %Std error of M9
SE_M10_dum = (Avar_hat_b_dum(11,11)/(n-1))^(1/2);      %Std error of M10
SE_M11_dum = (Avar_hat_b_dum(12,12)/(n-1))^(1/2);      %Std error of M11
SE_M12_dum = (Avar_hat_b_dum(13,13)/(n-1))^(1/2);      %Std error of M12

dum_beta_hat(2,1)
SE_R_dum = ((e_hat_dum-mean(e_hat_dum))'*(e_hat_dum-mean(e_hat_dum))/(n-13))^(1/2) %Std error Re

%ROBUST T-STAT
robust_t_beta1_OLS_dum = dum_beta_hat(1,1)/SE_tb1_dum

%Optional T-test
%   Process By: Daniela Osterrieder, Rutgers University, 2015; 
%   Only minor changes were made and commnents left intact
%H0: beta = 0 vs H1: beta>0 
critical_value = icdf('normal',0.95,0,1); 
if robust_t_beta1_OLS_dum<critical_value
    'based on a significance test with robust SEs, we fail to reject beta=0'
    
else
    'based on a significance test with robust SEs, we reject beta=0 in favour of the alternative hypothesis that beta>0'
end;
%   end Osterrieder, 2015

%Next, we add seasonal monthly dummies with a constant
X_w_dummiesC(1:223,1)=X(:,2);                          %create X matrix
X_w_dummiesC(1:223,2:13)=dummies;                      %add dummies
X_w_dummiesC(1:223,14)=ones(223,1);                    %add intercept

dumC_beta_hat = ((X_w_dummiesC'*X_w_dummiesC)^(-1))*(X_w_dummiesC'*y);

e_hat_dumC = y-X_w_dummiesC*dumC_beta_hat;                   %estimation of resids
R_sqr_dumC = 1-(e_hat_dumC'*e_hat_dumC)/((y-mean_y)'*(y-mean_y));   %fit using R^2

B_dumC = diag(e_hat_dumC.^2);                           %Create diag of res^2
S_hat_dumC = (X_w_dummiesC'*B_dumC*X_w_dummiesC)/(n-14);       %sum xixi'*e^2/n
Sxx_dumC = X_w_dummiesC'*X_w_dummiesC/(n-14);            %Variance matrix of X
Avar_hat_b_dumC = Sxx_dumC^(-1)*S_hat_dumC*Sxx_dumC^(-1);   %Asymtotic variance X

SE_tb1_dumC = (Avar_hat_b_dumC(1,1)/(n-1))^(1/2);        %Std error of tb1
SE_M1_dumC = (Avar_hat_b_dumC(2,2)/(n-1))^(1/2);         %Std error of M1
SE_M2_dumC = (Avar_hat_b_dumC(3,3)/(n-1))^(1/2);         %Std error of M2
SE_M3_dumC = (Avar_hat_b_dumC(4,4)/(n-1))^(1/2);         %Std error of M3
SE_M4_dumC = (Avar_hat_b_dumC(5,5)/(n-1))^(1/2);         %Std error of M4
SE_M5_dumC = (Avar_hat_b_dumC(6,6)/(n-1))^(1/2);         %Std error of M5
SE_M6_dumC = (Avar_hat_b_dumC(7,7)/(n-1))^(1/2);         %Std error of M6
SE_M7_dumC = (Avar_hat_b_dumC(8,8)/(n-1))^(1/2);         %Std error of M7
SE_M8_dumC = (Avar_hat_b_dumC(9,9)/(n-1))^(1/2);         %Std error of M8
SE_M9_dumC = (Avar_hat_b_dumC(10,10)/(n-1))^(1/2);       %Std error of M9
SE_M10_dumC = (Avar_hat_b_dumC(11,11)/(n-1))^(1/2);      %Std error of M10
SE_M11_dumC = (Avar_hat_b_dumC(12,12)/(n-1))^(1/2);      %Std error of M11
SE_M12_dumC = (Avar_hat_b_dumC(13,13)/(n-1))^(1/2);      %Std error of M12
SE_B0_dumC = (Avar_hat_b_dumC(14,14)/(n-1))^(1/2);      %Std error of intercept

SE_R_dumC = ((e_hat_dumC-mean(e_hat_dumC))'*(e_hat_dumC-mean(e_hat_dumC))/(n-13))^(1/2); %Std error Re

%NOTES
%We can see with the use of an intercept the Standard Errors explode
%Sxx is nearly singular with a recurring identity matrix and a vector of 1s

%Proof:
rcond(Sxx_dumC)

%The former dummy regression is so far the most significant in explaining the changes
%of our independent variable. The dummy with constant regression is, of
%course, not a reasonable regression judging by the standard errors and the
%dummy w/ constant trap addressed in other texts.

%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Part (j)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Start Part (k)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimate the Fama Regression for 1/53-7/71 using PA1 from the dataset
PAfama(:,2)=Mishkin(36:258,5);
PAfama(:,1)=ones;
ypafama=Mishkin(36:258,3);

PA1_fama_b_hat=((PAfama'*PAfama)^(-1))*(PAfama'*ypafama);
e_hat_fama = ypafama-PAfama*PA1_fama_b_hat;               %estimation of resids
mean_yfama=mean(ypafama);                                      %calc mean y
R_sqr_fama = 1-(e_hat_fama'*e_hat_fama)/((ypafama-mean_yfama)'*(ypafama-mean_yfama));   %R sqr for regression
n_fama=length(ypafama);
Bf = diag(e_hat_fama.^2);                                  %Create diag of res^2
S_hat_fama = (PAfama'*Bf*PAfama)/(n_fama-1);                               %sum xixi'*e^2/n
Sxx_fama = PAfama'*PAfama/(n_fama-1);                                     %Variance matrix of X
Avar_hat_betas_fama = Sxx_fama^(-1)*S_hat_fama*Sxx_fama^(-1);   %Asymtotic variance X

SE_B0_fama = (Avar_hat_betas_fama(1,1)/(n_fama-1))^(1/2);               %Std error of intcpt
SE_B1_fama = (Avar_hat_betas_fama(2,2)/(n_fama-1))^(1/2);               %Std error of Beta tb1
SE_R_fama = ((e_hat_fama-mean(e_hat_fama))'*(e_hat_fama-mean(e_hat_fama))/(n_fama-2))^(1/2) %Std error Re

%ROBUST T-STAT
robust_t_b1_OLS_fama = PA1_fama_b_hat(2,1)/SE_B1_fama

%Optional T-test
%   Process By: Daniela Osterrieder, Rutgers University, 2015; 
%   Only minor changes were made and commnents left intact
%H0: beta = 0 vs H1: beta>0 
critical_value = icdf('normal',0.95,0,1); 
if robust_t_b1_OLS_fama<critical_value
    'based on a significance test with robust SEs, we fail to reject beta=0'
    
else
    'based on a significance test with robust SEs, we reject beta=0 in favour of the alternative hypothesis that beta>0'
end;
%   end Osterrieder, 2015

%PLOTS & Tables
pa1_reg=figure;
set(pa1_reg,'name','Scatter and Best fit line FAMA REGRESSION PA1')
pa1_reg,
scatter(PAfama(:,2),ypafama);
hold on
plot(PAfama(:,2),PA1_fama_b_hat(2)*PAfama(:,2)+PA1_fama_b_hat(1));

R_tbl_fama(1,1:2)=PA1_fama_b_hat;
R_tbl_fama(2,1:2)=[SE_B0_fama,SE_B1_fama];
R_tbl_fama(3,3:6)=[R_sqr_fama,mean_yfama,SE_R_fama,robust_t_b1_OLS_fama];

freg_fama=figure('Position', [150 150 600 250]);
set(freg_fama, 'name', 'FAMA Regression Stats w/ PA1','numbertitle','off');

r2names_fama={'Coefficients','Std Error','Statistics'};
c2names_fama={,'B0','B1','R^2','Mean of y', 'SER in %','T-stat'};
Regression_table = uitable('Data',R_tbl_fama,...
    'RowName',r2names_fama,'ColumnName',c2names,'Tag',...
    'Fama Regression Inflation Rt by PA1, January 1953-July 1971',...
    'Parent', freg_fama,'Position',[40 40 550 170]);

%NOTES
%The notable differences are produced in the output. The difference in 
%coeficients and t statistic should be noted

%%%%%%%%%%%%%%%%%%%%%%%%%%%END PART (k)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Start Part (l)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimate the Regression for 10/79- using PA1 from the dataset
PAfama79(:,2)=Mishkin(358:end,5);
PAfama79(:,1)=ones;
ypafama79=Mishkin(358:end,3);

PA1_fama79_b_hat=((PAfama79'*PAfama79)^(-1))*(PAfama79'*ypafama79);
%!@#$%^&*()_+
e_hat_fama79 = ypafama79-PAfama79*PA1_fama79_b_hat;               %estimation of resids
mean_yfama79=mean(ypafama79);                                      %calc mean y
R_sqr_fama79 = 1-(e_hat_fama79'*e_hat_fama79)/((ypafama79-mean_yfama79)'*(ypafama79-mean_yfama79));   %R sqr for regression
n_fama79=length(ypafama79);
Bf79 = diag(e_hat_fama79.^2);                                  %Create diag of res^2
S_hat_fama79 = (PAfama79'*Bf79*PAfama79)/(n_fama79-1);                               %sum xixi'*e^2/n
Sxx_fama79 = PAfama79'*PAfama79/(n_fama79-1);                                     %Variance matrix of X
Avar_hat_betas_fama79 = Sxx_fama79^(-1)*S_hat_fama79*Sxx_fama79^(-1);   %Asymtotic variance X

SE_B0_fama79 = (Avar_hat_betas_fama79(1,1)/(n_fama79-1))^(1/2);               %Std error of intcpt
SE_B1_fama79 = (Avar_hat_betas_fama79(2,2)/(n_fama79-1))^(1/2);               %Std error of Beta tb1
SE_R_fama79 = ((e_hat_fama79-mean(e_hat_fama79))'*(e_hat_fama79-mean(e_hat_fama79))/(n_fama79-2))^(1/2) %Std error Re

%ROBUST T-STAT
robust_t_b1_OLS_fama79 = PA1_fama79_b_hat(2,1)/SE_B1_fama79

%Optional T-test
%   Process By: Daniela Osterrieder, Rutgers University, 2015; 
%   Only minor changes were made and commnents left intact
%H0: beta = 0 vs H1: beta>0 
critical_value = icdf('normal',0.95,0,1); 
if robust_t_b1_OLS_fama79<critical_value
    'based on a significance test with robust SEs, we fail to reject beta=0'
    
else
    'based on a significance test with robust SEs, we reject beta=0 in favour of the alternative hypothesis that beta>0'
end;
%   end Osterrieder, 2015

%PLOTS & Tables
pa1_reg79=figure;
set(pa1_reg79,'name','Scatter and Best fit line FAMA REGRESSION PA1 (Post Oct 79)')
pa1_reg,
scatter(PAfama79(:,2),ypafama79);
hold on
plot(PAfama79(:,2),PA1_fama79_b_hat(2)*PAfama79(:,2)+PA1_fama79_b_hat(1));

R_tbl_fama79(1,1:2)=PA1_fama79_b_hat;
R_tbl_fama79(2,1:2)=[SE_B0_fama79,SE_B1_fama79];
R_tbl_fama79(3,3:6)=[R_sqr_fama79,mean_yfama79,SE_R_fama79,robust_t_b1_OLS_fama79];

freg_fama79=figure('Position', [150 150 550 250]);
set(freg_fama79, 'name', 'FAMA Regression Stats w/ PA1 from OCT 79-','numbertitle','off');

r2names_fama79={'Coefficients','Std Error','Statistics'};
c2names_fama79={,'B0','B1','R^2','Mean of y', 'SER in %','T-stat'};
Regression_table = uitable('Data',R_tbl_fama79,...
    'RowName',r2names_fama,'ColumnName',c2names,'Tag',...
    'Fama Regression Inflation Rt by PA1, Oct 1979-',...
    'Parent', freg_fama79,'Position',[40 40 500 170]);

%NOTES
%The coefficient is lower, The y mean is much high and should be noted.
%Monetary policy of this time could explain outside influence to our model
%and our unusual coefficient.

%%%%%%%%%%%%%%%%%%%%%%%%%%%END PART (l)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%