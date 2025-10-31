# ivcrtest
Testing the validity of IVs using the Correlation Restriction

To install the package, use this code:
remotes::install_git("https://github.com/ratbekd/ivcrtest.git")
library(ivcrtest)

 One can test the package using this code:
 
 df <- wooldridge::card  # Use sample data set
 
data <- df[complete.cases(df), ]  # Remove missing values

attach(data)

H0 <- data.frame(lwage, educ,fatheduc, exper,expersq ,black ,south ,smsa ,smsa66,reg661, reg662 ,reg663, reg664,
                 reg665, reg666 ,reg667,reg668)
                 
Y <- as.character(colnames(H0))[1] ###OUTCOME variable

X <-as.character(colnames(H0))[2] ###ENDOEGENOUS variable

Z <-as.character(colnames(H0))[3] ###IV variable

H <-as.character(colnames(H0))[-(1:3)] ##### EXOGENOUS variables

result <- iv_cr_test(data,X,Y,H,Z,  n, k=-1,alpha = 0.05,seed = 123,rxu_range = c(0, 0.8),bias_mc = TRUE,mc_B = 500)

print(result)
