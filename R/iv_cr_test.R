#' Perform Instrumental Variable Correlation Restriction Test
#'
#' @param data A data frame containing all variables.
#' @param X Character. Name of the endogenous regressor.
#' @param Y Character. Name of the dependent variable.
#' @param H Character vector. Names of exogenous control variables.
#' @param Z Character. Name of the instrumental variable.
#' @param n Ignored; determined internally from \code{data}.
#' @param k Numeric scalar for direction of endogeneity. Default \code{-1}.
#' @param alpha Numeric significance level. Default \code{0.05}.
#' @param seed Integer seed for reproducibility. Default \code{123}.
#' @param rxu_range Numeric vector of length 2 for the \eqn{r_{xu}} grid.
#' @param bias_mc Logical. Monte Carlo bias correction. Default \code{FALSE}.
#' @param mc_B Integer number of MC draws. Default \code{500}.
#' @return A data frame with confidence intervals and probabilities.
#' @export
#' @importFrom stats complete.cases as.formula lm resid predict cor cov
#' @importFrom stats pf pchisq var df.residual qt sd rnorm pnorm qnorm
#' @importFrom stats optim quantile fitted
#' @examples
#' data <- wooldridge::card
#' data <- data[complete.cases(data), ]
#' Y <- "lwage"
#' X <- "educ"
#' Z <- "fatheduc"
#' H <- c("exper", "expersq", "black", "south", "smsa", "smsa66",
#'        "reg661", "reg662", "reg663", "reg664",
#'        "reg665", "reg666", "reg667", "reg668")
#' result <- iv_cr_test(data, X, Y, H, Z,
#'                      k = -1, alpha = 0.05, seed = 123,
#'                      rxu_range = c(0, 0.8))
#' print(result)
#'

iv_cr_test <- function(data,X,Y,H,Z,  n=NULL, k=-1,
                    alpha = 0.05,
                    seed = 123,
                    rxu_range = c(0, 0.8),
                    bias_mc = FALSE,
                    mc_B = 500) {

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
  ## ── keep only complete cases on variables we actually use ─────────────────
  vars_needed <- unique(c(Y, X, H, Z))
  data <- data[complete.cases(data[, vars_needed]), vars_needed]

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
  s <- min(length(y),length(x))
  x <- x[1:s]
  y <- y[1:s]
  xx <- min(length((x-mean(x))),length((y-mean(y))))
  data <- data[1:xx,]
  #saving the transformed x and y
  data$x<-(x-mean(x))
  data$y<-(y-mean(y))
  s <- min(length(y),length(x))
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


  df <- data.frame(x = x[1:k1], y = y[1:k1], z = z[1:k1])
  r_xz<-cor(df$x,df$z,use = "complete.obs",method=c("pearson"))
  r_yz<-cor(df$y,df$z,use = "complete.obs",method=c("pearson"))
  r_yx<-cor(df$y,df$x,use = "complete.obs",method=c("pearson"))
  # #

  ############################

  i=1
  cat("\n=== Running MCUB Method Membership Test ===\n")
  bias_target <- 0  # or a pre-computed bias value
  res <- check_compatibility(df, i=1, n = nrow(df), k = 1,
                             alpha = 0.05, rxu_range = c(0.0, 0.8))

  res

}
