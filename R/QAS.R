#'@title Quasi analytical solution for logit.
#'
#'@description \code{QAS.func} is used to replace numerical optimization with a quasi-analytical approach for logit models on big data. It returns the coefficients, predicted values and quality criteria for the provided variables.
#'
#'
#'@param frml an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#'@param data a data frame containing the variables in the model (or object coercible by \code{\link[base]{as.data.frame}} to a data frame).  Details of the structure of the data are given under 'Details'.
#'@param weight an optional vector of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. In case of NULL, each case is weighted with 1.
#'@param seed saving the state of a random process. Should be NULL or a numeric vector. In case of NULL a seed is generated at random.
#'
#'@details A typical predictor has the form dependent_Variable '~' independent_Variables.\cr The dependent_Variable has \strong{two} categories.\cr If there is more than one independent_Variable, they can be combined with a '+'.
#'@details The data frame can \strong{not} contain any missing values.\cr \strong{Metric} variables have to be of type \strong{numeric}. All \strong{other variables} have to be of type \strong{integer}.\cr The first variable in the dataset hat to be the dependent variable.\cr The scale of large numbers has to be reduced e.g. standardization.
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
QAS.func <- function(frml, data = NULL, weight = NULL, seed = NULL) {

   if (is.null(seed)) {
     seed <- floor(runif(1, 0, 1e6))
   }
   if (!is.null(seed)) {
    set.seed(seed)
   }

  # ----------------------------------------------------------- Datenaufbereitung

  model_data <- model.frame(frml, data = data)
  row.names(model_data) <- 1:nrow(data)

  if(is.integer(model_data[1])==FALSE) {

    model_data[1] <- as.integer(model_data[,1])

  }

  Y <- model_data[1]
  X <- model_data[-1]

  #---------------------------------- wenn keine Gewichtung gegeben, jeden Fall mit 1 gewichten

  if (is.null(weight)) {

    weight <- rep(1,nrow(model_data))
  }

  # --------------------------------  which Xvars have no variance

  if (sum(unlist(lapply(model_data, function(col) length(unique(col)))) == 1) == 1) {

    warning("remove variables with no variance.")
    out <- lapply(X, function(x) length(unique(x)))
    want <- which(!out > 1)
    X[,unlist(want)] <- NULL
    rm(out,want)

    out <- lapply(model_data, function(x) length(unique(x)))
    want <- which(!out > 1)
    model_data[,unlist(want)] <- NULL
    rm(out,want)
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
  catdat <- apply(numdat,2,categorize)
  catdatdat <- do.call(cbind,lapply(catdat,function(x)x[[1]]))
  model_data[,numindex] <- catdatdat
  means <- do.call(cbind,lapply(catdat,function(x)x[[2]]))[-1,]

  # -------------------------------- Cell definition

  # let's make 'Cell' in column going immediately after the data:
  model_data[,'Cell'] <- rep(0, nrow(model_data))

  # -------------------------------- Zeilen für contingency_function vorbereiten

  keeps <- names(X)
  model_data[,'Cell'] <- apply(as.matrix(model_data[keeps]), 1, function(x) paste(x, collapse = ' '))
  length(unique(model_data$Cell))
  rm(keeps)

  # ------------------------------- construct contingency y variable

  model_data <- contingency_function(data = model_data, model = Y)

  keeps <- c('LogOdds.PriorProb',names(X))
  modelling.all <- model_data[keeps]
  rm(keeps)

  # ------------------------------ Linear Regression

  #  ONE Steps to calculate:

  #  I.   calculate OLS.beta:
  #       OLS.beta = (X'X)^(-1)X'y

  results <- lm(formula = LogOdds.PriorProb ~. , data = modelling.all, na.action = na.exclude, weights = weight)

  # ----------------------------- Y_hat berechnen

  intercept <- rep(1,nrow(X))
  X <- cbind(intercept, X)

  #ist eingl kein X' !! sondern nur X
  X_strich <- (as.matrix(X))
  dim(X_strich)

  param <- t(as.matrix(results$coefficients))
  V <- tcrossprod(X_strich,param)
  eV <- exp(V)

  results$y_hat <- (eV/(1+eV))
  results$means_for_cat <- means

  results$seed <- seed

  results$fitted.values <- NULL
  results$residuals <- NULL
  results$assign <- NULL
  results$effects <- NULL
  results$rank <- NULL
  results$df.residual <- NULL
  results$xlevels <- NULL
  results$qr <- NULL

  # ---------------------------- AIC/BIC berechnen

  s <- (1/(nrow(data)-ncol(X)))*(sum((results$y_hat - Y)^2))

  results$AIC <- nrow(data)*log(s) + (2*ncol(X))

  results$BIC <- nrow(data)*log(s) + (ncol(X)*log(nrow(data)))

  return (results)

}

 #---------------------------- Wird nur intern verwendet, daher kein export
contingency_function <- function(data, model) {

  # --- --- ---
  Counts.of.Bought <- tapply(data[,1], data[,'Cell'], sum)
  Counts.of.Shown <- table(data[,'Cell'])

  CountBought <- as.numeric(Counts.of.Bought[match(data[,'Cell'], names(Counts.of.Bought))])
  CountsShown <- as.vector(Counts.of.Shown[match(data[,'Cell'], names(Counts.of.Shown))])
  PriorProb <- CountBought/CountsShown

  print(table(CountBought == 0))
  rm(Counts.of.Bought, Counts.of.Shown)
  rm(CountBought, CountsShown)

  # --- --- ---
  # change prob = 0 to eps, prob = 1 to (1 - eps):
  min.Prob <- min(PriorProb[PriorProb > 0])
  max.Prob <- max(PriorProb[PriorProb < 1])
  min.min.Prob <- min(min.Prob, 1 - max.Prob)
  eps <- min.min.Prob/100
  PriorProb[PriorProb == 0] <- eps
  PriorProb[PriorProb == 1] <- 1 - eps

  rm(min.Prob, max.Prob, min.min.Prob, eps)

  data$LogOdds.PriorProb <- log(PriorProb/(1 - PriorProb))

  return(data)
}

# -----------------------------------------------------
# K A T E G O R I S I E R U N G S - F U N K T I O N
#' @title Function for categorization after QAS.func
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
#'
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
#' @export
predictQAS <- function(QAS.res, data) {

  #QAS.res <- Result1
  #data <- data

  means_for_cut <- is.null(QAS.res$means_for_cat)

  if (means_for_cut == FALSE) {

    numvar <- names(which(sapply(data, class) == "numeric"))                 # alle numerics in data
    numvar <- numvar[numvar %in% c(attributes(QAS.res$terms)$term.labels)]   # alle die auch als UV werwendet wurden
    numdata <- as.data.frame(data[,numvar])
    colnames(numdata) <- numvar
    catdata <- matrix(NA, nrow = nrow(numdata), ncol = ncol(numdata))
    colnames(catdata) <- numvar

    if(length(numvar) == 1){
      for(j in numvar){

        QAS.res$means_for_cat <- sort(QAS.res$means_for_cat, decreasing = FALSE)
        catdata[,j] <- findInterval(numdata[,j],QAS.res$means_for_cat)       #Kategorisierung der Variablen

      }
    }

    if(length(numvar) > 1){
      for(j in numvar){

        catdata[,j] <- findInterval(numdata[,j],QAS.res$means_for_cat[,j])

      }

    }


    if(length(numvar) > 1){                                             # wenn mehr als eine Variable
      catdata <- do.call(cbind,lapply(catdata,function(x)x[[1]]))
    }

    data[,numvar] <- catdata                                            # Datenzusammenfügen
  }

  # y_hat berechnen

  Xnew1 <- c(attributes(QAS.res$terms)$term.labels)   # alle die auch als UV werwendet wurden
  Xnew <- as.data.frame(data[,Xnew1])
  intercept <- rep(1,nrow(Xnew))
  Xnew <- cbind(intercept, Xnew)

  #ist eingl kein X' !! sondern nur X
  X_strich <- (as.matrix(Xnew))
  param_QAS <- (as.matrix(QAS.res$coefficients))
  V <- (X_strich%*%param_QAS)
  eV <- exp(V)

  yhat_orig <- as.vector((eV/(1+eV)))

  return(list(yhat_orig=yhat_orig,data=data))
}
