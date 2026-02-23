# ivcrtest

# Summary of the CR test
We develop a Correlation Restriction (CR) test for assessing the validity of instrumental variables when the direction of endogeneity is known. Consistent with Gunsilius (2020), the test demonstrates that additional structural assumptions can render the instrument validity testable, even in models with continuous endogenous variables. Building on DiTraglia and Garcia-Jimeno (2021), our approach exploits the joint correlation structure among the instrument, regressor, and outcome to construct confidence intervals for partially identified parameters using a modified union-bound procedure. Monte Carlo simulations and empirical applications---including analyses of returns to education, recidivism, and development---demonstrate that the CR test reliably detects invalid instruments and characterises the range of exogeneity, thereby enabling formal empirical evaluation of instrument validity.

# The package
Testing the validity of IVs using the Correlation Restriction
The package computes the union bound CIs (simple and Bei (2024))  membership probability

Pr(r_zu=0 is a in C_n | r_xu in D), 

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

result <- iv_cr_test(data,X,Y,H,Z,  n, k=-1,alpha = 0.05,seed = 123,rxu_range = c(0, 0.8))

# Output

print( result)

One obtains these results:
                                 
Z                                   IV

plug_in                   [0.011, 0.806]

CI_Bei                    [-0.054, 0.845]

CI_simple                 [-0.108, 0.878]

Zero_in_CI_Bei                       ✓

p_zero                               1

n                                   1600
