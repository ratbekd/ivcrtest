#' Perform Instrumental Variable Correlation Restirction test
#' @param data Data frame
#' @param n Numeric scalar for size of the variables.
#' @param k Numeric scalar for direction of endogenity.
#' @param alpha Numeric scalar for the level of H0 test.
#' @param rxu_range Numeric interval for r_xu values.
#' @param seed Numeric scalar  for reproducibility.
#' @param bias_mc Logical value.
#' @param mc_B Numeric scalar for MC correction simulations.
#' @param data A data frame.
#' @param Y Name of the dependent variable.
#' @param X Name of the endogenous variable.
#' @param Z Name of the instrumental variable.
#' @param H Vector of exogenous variable names.
#' @return A list containing CIs and coverage or containment probabilities.
#' @export
#' @importFrom stats  predict
#' @importFrom stats  complete.cases
#' @importFrom stats  fitted
#' @importFrom stats pf
#' @importFrom stats pchisq
#' @importFrom stats var
#' @importFrom stats as.formula
#' @importFrom stats lm
#' @importFrom stats cor
#' @importFrom stats cov
#' @importFrom stats residuals
#' @importFrom stats df.residual
#' @importFrom stats qt
#' @importFrom stats sd
#' @importFrom stats rnorm
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats predict
#' @importFrom stats df
#' @importFrom stats optim
#' @importFrom stats resid
#' @importFrom stats quantile
#'
#' @examples
#' df <- wooldridge::card  # Use sample data set
#' data <- df[complete.cases(df), ]  # Remove missing values
#' attach(data)
#' H0 <- data.frame(lwage, educ,fatheduc, exper,expersq ,black ,south ,smsa ,smsa66,reg661, reg662 ,reg663, reg664,
#'reg665, reg666 ,reg667,reg668)
# Vector of variables used in the regression. Firts var is outcome, the second var is endogenous regressor
#' Y <- as.character(colnames(H0))[1] ###OUTCOME variable
#' X <-as.character(colnames(H0))[2] ###ENDOEGENOUS variable
#' Z <-as.character(colnames(H0))[3] ###IV variable
#' H <-as.character(colnames(H0))[-(1:3)] ##### EXOGENOUS variables
<<<<<<< HEAD
#' result <- iv_cr_test(data,X,Y,H,Z,  n, k=-1,alpha = 0.05,seed = 123,rxu_range = c(0, 0.8))
#' print(result)
=======
#' result <- iv_cr_test(data,X,Y,H,Z,  n, k=-1,alpha = 0.05,seed = 123,rxu_range = c(0, 0.8),bias_mc = TRUE,mc_B = 500)
#' results_list <- result
#' summary_df <- do.call(rbind, results_list)
# Force it to be a data frame
#' summary_df <- as.data.frame(summary_df)
#' summary_df
>>>>>>> 102ddd71075d4e85a96573aacee4e0717e4e90b5

iv_cr_test <- function(data,X,Y,H,Z,  n=NULL, k=-1,
                    alpha = 0.05,
                    seed = 123,
                    rxu_range = c(0, 0.8),
                    bias_mc = FALSE,
                    mc_B = 500) {
  library(dplyr)
  library(ivreg)
  library(lmtest)
  library(sandwich)
  library(AER)
  library(MASS)
  library(e1071)
  library(zoo)
  ####################################
  if (!all(c(Y, X, H,Z) %in% names(data))) {
    stop("Some variables are missing in the dataset.")
  }
  # Y <- as.character(colnames(H0))[1] ###OUTCOME variable
  # X <-as.character(colnames(H0))[2] ###ENDOEGENOUS variable
  # Z <-as.character(colnames(H0))[3] ###IV variable
  # H<- as.character(colnames(H0))[-(1:3)] ##### EXOGENOUS variables
  #print(H)
  # formula_str <- paste(paste0(Y," ~ ", X,"+"), paste(H, collapse = " + "))## Construct formula as a string
  # formula <- as.formula(formula_str)# Convert to formula object
  # #print(formula)

  # determine the number of rows
  n<-nrow(data)

  ### Initialization of the variables
  fitx=0
  fity=0
  fitz=0
  # ### the outcome variable
  ## Factoring out the effects of other exogenous variables
  formula_str <- paste(paste0(Y," ~ "), paste(H, collapse = " + "))## Construct formula as a string
  formula <- as.formula(formula_str)# Convert to formula object
  fity<-lm(formula, data=data)
  y<-resid(fity)
  # the endogenous variable
  ## Factoring out the effects of other exogenous variables
  formula_str <- paste(paste0(X," ~ "), paste(H, collapse = " + "))## Construct formula as a string
  formula <- as.formula(formula_str)# Convert to formula object
  fitx<-lm(formula, data=data)
  x<-resid(fitx)
  s <- max(length(y),length(x))
  x <- x[1:s]
  y <- y[1:s]
  xx <- min(length((x-mean(x))),length((y-mean(y))))
  data <- data[1:xx,]
  #saving the transformed x and y
  data$x<-(x-mean(x))
  data$y<-(y-mean(y))
  s <- max(length(y),length(x))
  x <- x[1:s]
  y <- y[1:s]
  options(digits=5)
  ### IV ######################################

  #z<-instrument # to use cor(xz)>0
  formula_str <- paste(paste0(Z," ~ "), paste(H, collapse = " + "))## Construct formula as a string
  formula <- as.formula(formula_str)# Convert to formula object
  fitz<-lm(formula, data=data)
  z<-resid(fitz)
  k1<-min(length(z),length(x),length(y))

  z<-predict(lm(z[1:k1]~x[1:k1]+y[1:k1],data=data))
  k1<-min(length(z),length(x),length(y))

  z<-predict(lm(z[1:k1]~x[1:k1]+y[1:k1],data=data))

  #df1 <- data.frame(x[1:k1],y[1:k1],z0p[1:k1])
  df <- data.frame(x = x[1:k1], y = y[1:k1], z = z[1:k1])
  r_xz<-cor(df$x,df$z,use = "complete.obs",method=c("pearson"))
  r_yz<-cor(df$y,df$z,use = "complete.obs",method=c("pearson"))
  r_yx<-cor(df$y,df$x,use = "complete.obs",method=c("pearson"))
  # #

  ############################

  i=1
  cat("\n=== Running MCUB Method Membership Test ===\n")
  bias_target <- 0  # or a pre-computed bias value
  res <- check_compatibility(df, 1, n = nrow(df), k = 1,
                             alpha = 0.05, rxu_range = c(0.0, 0.8))

  res

}
