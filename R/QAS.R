#'@title Quasi analytical solution for logit.
#'
#'@description \code{QAS.func} is used to replace numerical optimization with a quasi-analytical approach for logit models on big data. It returns the coefficients, predicted values and quality criteria for the provided variables.
#'
#'
#'@param frml an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.
#'@param data a data frame containing the variables in the model (or object coercible by \code{\link[base]{as.data.frame}} to a data frame).  Details of the structure of the data are given under 'Details'.
#'@param weights an optional vector of prior weights to be used in the fitting process. Should be NULL or a numeric vector. In case of NULL, each case is weighted with 1.
#'@param seed saving the state of a random process. Should be NULL or a numeric vector. In case of NULL a seed is generated at random.
#'@param tau an optional parameter proposed by King and Zeng (2001) which comprises prior information about the fraction of ones in the population of the dependent variable. It has to lie between 0 and 1.
#'
#'@details A typical predictor has the form dependent_Variable '~' independent_Variables.\cr The dependent_Variable has \strong{two} categories.\cr If there is more than one independent_Variable, they can be combined with a '+'.
#'@details The data frame must \strong{not} contain any missing values.\cr \strong{Metric} variables have to be of type \strong{numeric}. All \strong{other variables} have to be of type \strong{integer}.\cr The first variable in the dataset hat to be the dependent variable.\cr The scale of large numbers has to be reduced e.g. standardization.
#'
#'@return An object of class \emph{QAS.func} is a list containing at least the following components:
#'
#'\describe{
#'  \item{\code{coefficients}}{a vector of coefficients}
#'  \item{\code{weights}}{the working weights}
#'  \item{\code{call}}{the call of the final function within QAS.func}
#'  \item{\code{terms}}{the term object used}
#'  \item{\code{model}}{the model frame}
#'  \item{\code{means.for.cat}}{the cut points of the metric variables for a categorization of the original dataset}
#'  \item{\code{categorized.variables}}{the variables that have been categorized within QAS}
#'  \item{\code{seed}}{used seed for calculations}
#'
#'}
#'
#'@author This method is based on the research work of Stan Lipovetsky and Birgit Stoltenberg.
#'
#'@references
#'King, G. & Zeng, L. (2001), Logistic Regression in Rare Events Data, \emph{Political Analysis}, No. 9 / 2001 \cr
#'Lipovetsky, S. (2014), Analytical closed-form solution for binary logit regression by categorical predictors, \emph{Journal of Applied Statistics}, No. 42 / 2015 \cr
#'Lipovetsky, S. & Conklin, M. (2014), Best-Worst Scaling in analytical closed-form solution, \emph{The Journal of Choice Modelling}, No. 10 / 2014 \cr
#'Stoltenberg, B. (2016), Using logit on big data - from iterative methods to analytical solutions, \emph{GfK Verein Working Paper Series}, No. 3 / 2016 \cr
#'
#'@examples
#' # generate Data
#' y <- as.integer(c(1,0,0,0,1,1,1,0,0,1))
#' x <- c(15,88,90,60,24,30,26,57,69,18)
#' z <- as.integer(c(3,2,2,1,3,3,2,1,1,3))
#' example_data <- data.frame(y,x,z)
#'
#' # deploy QAS.func-Function
#' result <- QAS.func(y~x+z, data=example_data, weights=NULL, seed=NULL, tau = NULL)
#'
#'@importFrom stats runif lm model.frame
#'
#'@export
QAS.func <- function(frml, data = data, weights = NULL, seed = NULL, tau = NULL) {
  NAtest <-  unlist(lapply(data, anyNA))

  if (any(NAtest)){
    NAnames <- names(NAtest[NAtest == TRUE])
    stop(paste("These variables comprise missing values:",
               paste(NAnames, collapse = ", ")))
    rm(NAnames)
  }
  rm(NAtest)

  if (is.null(seed)) {
    seed <- floor(runif(1, 0, 1e6))
  } else {
    set.seed(seed)
  }

  # ----------------------------------------------------------- Datenaufbereitung

  model_data <- model.frame(frml, data = data)

  if(is.integer(model_data[,1])==FALSE) {

    stop("First variable is no integer.")
  }

  z <- rle(sort(model_data[,1]))
  if(length(z[["values"]])!=2) {

    stop("Dependent variable does not have two categories.")
  }

  fact <- vapply(model_data,is.factor,c(is.factor=FALSE))

  if (sum(fact)>0) {

    print(names(Filter(is.factor, model_data)))
    stop("Factors in the dataset")
  }

  char <- vapply(model_data,is.character,c(is.character=FALSE))

  if (sum(char)>0) {

    print(names(Filter(is.character, model_data)))
    stop("Characters in the dataset")
  }

  if (!is.null(tau)) {

    if(tau <= 0 | tau >= 1) {
      stop("Tau has to lie between 0 and 1.")
    }
  }

  nr <- nrow(model_data)
  nc <- ncol(model_data)

  # --------------------------------  which X do not have any variance

  want <- apply(model_data,2,function(x){diff(range(x))}) >0
  swant <- sum(want)
  if(swant < nc){

    model_data <- model_data[,want]
    warning(paste(nc-swant,"variable(s) with no variance removed."))
  }

  # -------------------------------- Kategorisierung metrischer Variablen

  categorize <- function(y){
    my <- mean(y)
    my2 <- tapply(y,y > my,mean)
    rulez <- c(-.Machine$double.xmax,my2[1],my,my2[2])
    rulez <- sort(rulez, decreasing = FALSE)
    caty <- findInterval(y,rulez)
    return(list(caty,rulez))
  }

  numindex <- which(sapply(model_data, function(x)length(unique(x))) > 6)
  
  if (length(numindex) != 0) {

  numdat <- as.data.frame(model_data[,numindex])

  colnames(numdat) <- names(numindex)
  catdat <- apply(numdat,2,categorize)
  catdatdat <- do.call(cbind,lapply(catdat,function(x)x[[1]]))
  model_data[,numindex] <- catdatdat

  if (length(numindex) == 1){
    model_data[,numindex] <- as.integer(model_data[,numindex])
  }

  means <- do.call(cbind,lapply(catdat,function(x)x[[2]]))[-1,,drop=FALSE]
  rownames(means) <- c(1:nrow(means))
  
  }

  # -------------------------------- Cell definition

  Y <- model_data[1]
  X <- model_data[-1]

  # let's make 'Cell' in column going immediately after the data:
  model_data[,'Cell'] <- rep(0, nr)

  # -------------------------------- Zeilen fuer contingency_function vorbereiten

  keeps <- names(X)
  model_data[,'Cell'] <- apply(as.matrix(model_data[keeps]), 1, function(x) paste(x, collapse = ' '))
  rm(keeps)

  Cell <- model_data[,'Cell']

  # ------------------------------- construct contingency y variable

  LogOdds.PriorProb <- contingency_function(Y = Y, Cell = Cell)

  model_data <- cbind(LogOdds.PriorProb,model_data)
  keeps <- c('LogOdds.PriorProb',names(X))
  modelling.all <- model_data[keeps]
  rm(keeps)

  # ------------------------------ Linear Regression

  #  Calculate OLS.beta:
  #       OLS.beta = (X'X)^(-1)X'y


  if(is.null(weights)) {

    results <- lm(formula = LogOdds.PriorProb ~. , data = modelling.all)

  } else {

    results <- lm(formula = LogOdds.PriorProb ~. , data = modelling.all, weights = weights)
  }

  # NA-Parameter durch "0" ersetzen

  nacoef <- which(is.na(results$coefficients))
  results$coefficients[nacoef] <- 0

  # Adapt intercept

  if (!is.null(tau)) {
    ymean <- mean(Y[,1])
    results$coefficients[1] <- (results$coefficients[1])-log(((1-tau)/tau)*(ymean/(1-ymean)))
    rm(ymean)
  }

  if (length(numindex) != 0) {
  results$means.for.cat <- means
  results$categorized.variables <- colnames(means)
  }
  results$seed <- seed
  # apply original y to model output
  results$model <- cbind(data[1], results$model[-1])

  results$residuals <- NULL
  results$fitted.values <- NULL
  results$assign <- NULL
  results$effects <- NULL
  results$rank <- NULL
  results$df.residual <- NULL
  results$xlevels <- NULL
  results$qr <- NULL

  return (results)

}

#---------------------------- Wird nur intern verwendet, daher kein Export
contingency_function <- function(Y, Cell) {

  Counts.of.Bought <- tapply(Y[,1], Cell, mean)

  PriorProb <- as.numeric(Counts.of.Bought[match(Cell, names(Counts.of.Bought))])

  rm(Counts.of.Bought)

  # --- --- ---
  # change prob = 0 to eps, prob = 1 to (1 - eps):
  min.Prob <- min(PriorProb[PriorProb > 0])
  max.Prob <- max(PriorProb[PriorProb < 1])
  min.min.Prob <- min(min.Prob, 1 - max.Prob)
  eps <- min.min.Prob/100
  PriorProb[PriorProb == 0] <- eps
  PriorProb[PriorProb == 1] <- 1 - eps

  rm(min.Prob, max.Prob, min.min.Prob, eps)

  LogOdds.PriorProb <- log(PriorProb/(1 - PriorProb))

  return(LogOdds.PriorProb)
}

# -----------------------------------------------------
# K A T E G O R I S I E R U N G S - F U N K T I O N
#'@title Function for categorization after QAS.func
#'
#'@description \code{predictQAS} is used to categorize  numeric variables in a dataset the same way as the QAS.func-Function for calculating the coefficients.
#'
#'
#'@param QAS.res The results of the QAS.func-Function
#'@param data A data frame containing the variables in the model for QAS.func
#'@param same a logical vector which indicates whether prediction of the dependent variable is applied for the same dataset used for QAS.func. Its default is FALSE.
#'
#'@return An object of class \emph{predictQAS} is a list containing at least the following components:
#'
#'\describe{
#'  \item{\code{yhat}}{the predicted values of the dependent variable based on the coefficients of QAS.func and the dataset used for categorization}
#'  \item{\code{transformed_data}}{the dataset with categorized numeric variables}
#'  \item{\code{AIC}}{A version of Akaike's Information Criterion}
#'  \item{\code{BIC}}{A version if Bayesian Information Criterion}
#'  }
#'
#'@examples
#' # generate Data
#' y <- as.integer(c(1,0,0,0,1,1,1,0,0,1))
#' x <- c(15,88,90,60,24,30,26,57,69,18)
#' z <- as.integer(c(3,2,2,1,3,3,2,1,1,3))
#' example_data <- data.frame(y,x,z)
#'
#' # deploy QAS.func-Function
#' result1 <- QAS.func(y~x+z, data=example_data, weights=NULL, seed=NULL, tau = NULL)
#'
#' # deploy predictQas-Function
#' result2  <- predictQAS(QAS.res = result1, data=example_data)
#'
#'@export
predictQAS <- function(QAS.res, data, same = FALSE) {

  cat.var <- QAS.res$categorized.variables

  if (same){
    transformed_data <- QAS.res$model
  } else {
    transformed_data <- data

    if (!is.null(cat.var)) {

      cat.data <- data[names(data) %in% cat.var]
      means <- QAS.res$means.for.cat

      for(j in cat.var){

        cat.data[,j] <- findInterval(cat.data[,j],means[,j])+1 #for same categorization as in QAS.func
      }
      transformed_data[names(transformed_data) %in% cat.var] <- cat.data
    }
  }

  # ----------------------------- Y_hat berechnen

  Y <- transformed_data[1]
  X <- transformed_data[-1]
  intercept <- rep(1,nrow(X))
  X <- as.matrix(cbind(intercept, X))

  param <- t(as.matrix(QAS.res$coefficients))
  V <- tcrossprod(X,param)
  eV <- exp(V)

  yhat <- as.vector(eV/(1+eV))

  # ---------------------------- AIC/BIC berechnen

  nr <- nrow(transformed_data)
  nc <- ncol(transformed_data)

  logL <- sum(Y*log(yhat)) + sum((1-Y)*log(1-yhat))

  AIC <- 2*nr - 2*logL

  BIC <- log(nc)*nr - 2*logL

  return(list(yhat=yhat,transformed_data=transformed_data, AIC = AIC,
              BIC = BIC))

}

# -----------------------------------------------------
# QAS for MICE - F U N K T I O N

#'@title Function to use QAS in Mice
#'
#'@description \code{mice.impute.QAS} is used to prepare the \code{predictQAS} output for mice
#'
#'@param y The target variable
#'@param ry a logical vector, which indicates the missing values in the target variable
#'@param x a data frame with the independent variables
#'@param boot if TRUE, a bayesian bootstrap sample is applied
#'@param ... further objects from other functions
#'
#'
#'@return An object of class \emph{mice.impute.QAS} is a vector containing the imputed values of the target variable.
#'
#'@examples
#' # generate Data
#' y <- as.integer(c(1,0,0,0,1,1,1,0,0,1))
#' x <- c(15,88,90,60,24,30,26,57,69,18)
#' z <- as.integer(c(3,2,2,1,3,3,2,1,1,3))
#' example_data <- data.frame(y,x,z)
#'
#' # deploy QAS.func-Function
#' result1 <- QAS.func(y~x+z, data=example_data, weights=NULL, seed=NULL, tau=NULL)
#'
#' # deploy predictQas-Function
#' result2  <- predictQAS(QAS.res = result1, data=example_data)
#'
#' # generate logical vector and data frame with independent variables
#' ry <- c(TRUE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE)
#' x <- data.frame(x,z)
#'
#' # deploy mice.impute.QAS-Funktion
#' impute <- mice.impute.QAS(y=y, ry=ry, x=example_data, boot=TRUE)
#'
#' # run mice function
#' library(mice)
#' final <- mice(nhanes2,method=c("","pmm","QAS","pmm"))
#'
#'@importFrom stats formula rbinom
#'
#'@export

mice.impute.QAS <- function (y, ry, x, boot=TRUE, ...) {
  ty <- names(table(y))
  if(length(ty)==2) {
    y.QAS <- as.integer(y==ty[2])
  } else if(length(ty)==1){
    y.QAS <- as.integer(y==ty[1])
  } else {
    stop("target variable has more than two categories")
  }

  data <- data.frame(y.QAS,x)
  colnames(data) <- paste0("v",c(1:(ncol(x)+1)))

  nobs <- sum(ry)
  nmis <- sum(!ry)

  data.obs <- data[ry, ]
  y.obs <- data.obs[,1]

  frml <- formula(paste("v1~",paste(paste0("v",c(2:(ncol(x)+1))),collapse="+")))


  if (diff(range(y.obs)) == 0){

    yimp <- rep(y.obs[1],nmis)
  } else {
    if(boot){
      random <- runif(nobs - 1)
      sorted <- c(0, sort(random), 1)
      weights <- diff(sorted) * nobs
      rm(random, sorted)

      QAS.res <- QAS.func(frml = frml,data = data.obs, weights = weights)
    } else {
      QAS.res <- QAS.func(frml = frml,data = data.obs)
    }

    yhat.mis <- predictQAS(QAS.res, data[!ry,])

    yimp <- sapply(yhat.mis$yhat,rbinom,size=1,n=1)

  }

  yimp[yimp==1] <- ty[2]
  yimp[yimp==0] <- ty[1]

  if (class(y) == "factor"){
    yimp <- as.factor(yimp)
  } else {
    class(yimp) <- class(y)
  }

  return(yimp)
}
