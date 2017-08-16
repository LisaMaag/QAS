#'@title Quasi analytical solution for logit.
#'
#'@description \code{QAS.func} is used to replace numerical optimization with a quasi-analytical approach for logit models on big data. It returns the coefficients, predicted values and quality criteria for the provided variables.
#'
#'
#'@param frml an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#'@param data a data frame containing the variables in the model (or object coercible by \code{\link[base]{as.data.frame}} to a data frame).  Details of the structure of the data are given under 'Details'.
#'@param weights an optional vector of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. In case of NULL, each case is weighted with 1.
#'@param seed saving the state of a random process. Should be NULL or a numeric vector. In case of NULL a seed is generated at random.
#'
#'@details A typical predictor has the form dependent_Variable '~' independent_Variables.\cr The dependent_Variable has \strong{two} categories.\cr If there is more than one independent_Variable, they can be combined with a '+'.
#'@details The data frame must \strong{not} contain any missing values.\cr \strong{Metric} variables have to be of type \strong{numeric}. All \strong{other variables} have to be of type \strong{integer}.\cr The first variable in the dataset hat to be the dependent variable.\cr The scale of large numbers has to be reduced e.g. standardization.
#'
#'
#'@return An object of class \emph{QAS.func} is a list containing at least the following components:
#'
#'\describe{
#'  \item{\code{coefficients}}{a vector of coefficients}
#'  \item{\code{weights}}{the working weights}
#'  \item{\code{call}}{the call of the final function within QAS.func}
#'  \item{\code{terms}}{the term object used}
#'  \item{\code{model}}{the model frame}
#'  \item{\code{y_hat}}{the predicted y values}
#'  \item{\code{means_for_cat}}{the cut points of the metric variables for a categorization of the original dataset}
#'  \item{\code{seed}}{used seed for calculations}
#'  \item{\code{AIC}}{A version of Akaike's Information Criterion}
#'  \item{\code{BIC}}{A version if Bayesian Information Criterion}
#'}
#'
#'@author This method is based on the research work of Stan Lipovetsky and Birgit Stoltenberg.
#'
#'@references
#'Lipovetsky, S. (2014), Analytical closed-form solution for binary logit regression by categorical predictors, \emph{Journal of Applied Statistics}, No. 42 / 2015 \cr
#'Lipovetsky, S. & Conklin, M. (2014), Best-Worst Scaling in analytical closed-form solution, \emph{The Journal of Choice Modelling}, No. 10 / 2014 \cr
#'Stoltenberg, B. (2016), Using logit on big data – from iterative methods to analytical solutions, \emph{GfK Verein Working Paper Series}, No. 3 / 2016 \cr
#'
#'@examples
#' # generate Data
#' y <- as.integer(c(1,0,0,0,1,1,1,0,0,1))
#' x <- c(15,88,90,60,24,30,26,57,69,18)
#' z <- as.integer(c(3,2,2,1,3,3,2,1,1,3))
#' example_data <- data.frame(y,x,z)
#'
#' # deploy QAS.func-Function
#' result <- QAS.func(y~x+z, data=example_data,weight=NULL,seed=NULL)
#'
#'
#'@importFrom stats runif lm model.frame na.exclude
#'
#'@export
QAS.func <- function(frml, data = data, weights = NULL, seed = NULL) {

  if (is.null(seed)) {
    seed <- floor(runif(1, 0, 1e6))
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ----------------------------------------------------------- Datenaufbereitung

  model_data <- model.frame(frml, data = data)

  if(is.integer(model_data[,1])==FALSE) {

    stop("first variable is no integer")
  }

  z <- rle(sort(model_data[,1]))
  if(length(z[["values"]])!=2) {

    stop("dependent variable does not have two categories")
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

  nr <- nrow(model_data)
  nc <- ncol(model_data)

  # --------------------------------  which Xvars have no variance

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

  numindex <- which(sapply(model_data, class) == "numeric")
  numdat <- as.data.frame(model_data[,numindex])
  colnames(numdat) <- names(numindex)
  catdat <- apply(numdat,2,categorize)
  catdatdat <- do.call(cbind,lapply(catdat,function(x)x[[1]]))
  model_data[,numindex] <- catdatdat
  means <- do.call(cbind,lapply(catdat,function(x)x[[2]]))[-1,,drop=FALSE]

  # -------------------------------- Cell definition

  Y <- model_data[1]
  X <- model_data[-1]

  # let's make 'Cell' in column going immediately after the data:
  model_data[,'Cell'] <- rep(0, nr)

  # -------------------------------- Zeilen für contingency_function vorbereiten

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

  # ------------------------------ Lineare Regression

  #  ONE Steps to calculate:
  #  I.   calculate OLS.beta:
  #       OLS.beta = (X'X)^(-1)X'y


  if(is.null(weights)) {

    results <- lm(formula = LogOdds.PriorProb ~. , data = modelling.all, na.action = na.exclude)

  } else {

    results <- lm(formula = LogOdds.PriorProb ~. , data = modelling.all, na.action = na.exclude, weights = weights)

  }

  # NA-Parameter durch "0" ersetzen

  nacoef <- which(is.na(results$coefficients))
  results$coefficients[nacoef] <- 0

  # ----------------------------- Y_hat berechnen

  intercept <- rep(1,nrow(X))
  X <- cbind(intercept, X)

  #ist eingl kein X' !! sondern nur X
  X_strich <- (as.matrix(X))

  param <- t(as.matrix(results$coefficients))
  V <- tcrossprod(X_strich,param)
  eV <- exp(V)

  results$fitted.values <- (eV/(1+eV))
  results$means_for_cat <- means
  results$seed <- seed

  results$residuals <- NULL
  results$assign <- NULL
  results$effects <- NULL
  results$rank <- NULL
  results$df.residual <- NULL
  results$xlevels <- NULL
  results$qr <- NULL

  # ---------------------------- AIC/BIC berechnen

  logL <- sum(Y*log(results$fitted.values)) + sum((1-Y)*log(1-results$fitted.values))

  results$AIC <- 2*nr - 2*logL

  results$BIC <- log(nc)*nr - 2*logL

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
#'
#'@return An object of class \emph{predictQAS} is a list containing at least the following components:
#'
#'\describe{
#'  \item{\code{yhat_orig}}{the predicted values of the dependent variable based on the coefficients of QAS.func and the dataset used for categorization}
#'  \item{\code{data}}{the dataset with categorized numeric variables}
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
#' result1 <- QAS.func(y~x+z, data=example_data,weight=NULL,seed=NULL)
#'
#' # deploy predictQas-Function
#' result2  <- predictQAS(QAS.res = result1, data=example_data)
#'
#'
#'@export
predictQAS <- function(QAS.res, data) {

  means_for_cut <- is.null(QAS.res$means_for_cat)

  if (means_for_cut == FALSE) {

    numvar <- names(which(sapply(data, class) == "numeric"))                 # alle numerics in data
    numvar <- numvar[numvar %in% c(attributes(QAS.res$terms)$term.labels)]   # alle die auch als UV verwendet wurden
    numdata <- as.data.frame(data[,numvar])
    colnames(numdata) <- numvar
    catdata <- matrix(NA, nrow = nrow(numdata), ncol = ncol(numdata))
    colnames(catdata) <- numvar


    for(j in numvar){

      QAS.res$means_for_cat[,j] <- sort(QAS.res$means_for_cat[,j], decreasing = FALSE)
      # QAS.res$means_for_cat <- sort(QAS.res$means_for_cat, decreasing = FALSE)

      catdata[,j] <- findInterval(numdata[,j],QAS.res$means_for_cat[,j])
      # catdata <- findInterval(numdata,QAS.res$means_for_cat)

    }

    data[,numvar] <- catdata                                            # Datenzusammenfügen
  }

  # y_hat berechnen

  Xnew1 <- c(attributes(QAS.res$terms)$term.labels)   # alle die auch als UV verwendet wurden
  Xnew <- as.data.frame(data[,Xnew1])
  intercept <- rep(1,nrow(Xnew))
  Xnew <- cbind(intercept, Xnew)

  #ist eingl kein X' !! sondern nur X
  X_strich <- (as.matrix(Xnew))
  param_QAS <- t(as.matrix(QAS.res$coefficients))
  V <- tcrossprod(X_strich,param_QAS)
  eV <- exp(V)

  yhat_orig <- as.vector((eV/(1+eV)))
  transformed_data <- data
  row.names(transformed_data) <- 1:nrow(transformed_data)            # Nummerierung der Fälle

  return(list(yhat_orig=yhat_orig,transformed_data=transformed_data))
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
#'@param boot if TRUE, a bootstrap sample ?
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
#' result1 <- QAS.func(y~x+z, data=example_data,weight=NULL,seed=NULL)
#'
#' # deploy predictQas-Function
#' result2  <- predictQAS(QAS.res = result1, data=example_data)
#'
#' # generate logical vector and data frame with independent variables
#' ry <- c(TRUE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE)
#' x <- data.frame(x,z)
#'
#' # deploy mice.impute.QAS-Funktion
#' impute <- mice.impute.QAS(y=y,ry=ry,x=example_data,boot=TRUE)
#'
#' # run mice function
#' library(mice)
#' final <- mice(nhanes2,method=c("","pmm","QAS","pmm"))
#'
#'@importFrom stats formula rbinom
#'
#'@export

mice.impute.QAS <- function (y, ry, x,boot=TRUE, ...) {

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
      weight <- diff(sorted) * nobs
      rm(random, sorted)

      QAS.res <- QAS.func(frml = frml,data = data.obs, weights = weight)
    } else {
      QAS.res <- QAS.func(frml = frml,data = data.obs)
    }

    yhat.mis <- predictQAS(QAS.res, data[!ry,])

    yimp <- sapply(yhat.mis$yhat_orig,rbinom,size=1,n=1)

  }

  yimp[yimp==1] <- ty[2]
  yimp[yimp==0] <- ty[1]

  class(yimp) <- class(y)

  return(yimp)
}
