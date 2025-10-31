# ivcrtest

# Summary of the CR test
We develop a Correlation Restriction (CR) test for assessing the validity of instrumental variables when the direction of endogeneity is known. Consistent with \citet{gunsilius2020}, the test shows that additional structural assumptions can render instrument validity testable even in models with continuous endogenous variables. Building on \citet{DiTraglia2021}, our approach exploits the joint correlation structure among the instrument, regressor, and outcome to construct confidence intervals for partially identified parameters using a modified union-bound procedure. Monte Carlo simulations and empirical applications---including analyses of returns to education, recidivism, and development---demonstrate that the CR test reliably detects invalid instruments and characterizes the range of exogeneity, thereby enabling formal empirical evaluation of instrument validity.

# The package
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
