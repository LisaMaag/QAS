% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QAS.R
\name{QAS.func}
\alias{QAS.func}
\title{Quasi analytical solution for logit.}
\usage{
QAS.func(frml, data = data, weights = NULL, seed = NULL, tau = NULL)
}
\arguments{
\item{frml}{an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.}

\item{data}{a data frame containing the variables in the model (or object coercible by \code{\link[base]{as.data.frame}} to a data frame).  Details of the structure of the data are given under 'Details'.}

\item{weights}{an optional vector of prior weights to be used in the fitting process. Should be NULL or a numeric vector. In case of NULL, each case is weighted with 1.}

\item{seed}{saving the state of a random process. Should be NULL or a numeric vector. In case of NULL a seed is generated at random.}

\item{tau}{an optional parameter proposed by King and Zeng (2001) which comprises prior information about the fraction of ones in the population of the dependent variable. It has to lie between 0 and 1.}
}
\value{
An object of class \emph{QAS.func} is a list containing at least the following components:

\describe{
 \item{\code{coefficients}}{a vector of coefficients}
 \item{\code{weights}}{the working weights}
 \item{\code{call}}{the call of the final function within QAS.func}
 \item{\code{terms}}{the term object used}
 \item{\code{model}}{the model frame}
 \item{\code{means.for.cat}}{the cut points of the metric variables for a categorization of the original dataset}
 \item{\code{categorized.variables}}{the variables that have been categorized within QAS}
 \item{\code{seed}}{used seed for calculations}

}
}
\description{
\code{QAS.func} is used to replace numerical optimization with a quasi-analytical approach for logit models on big data. It returns the coefficients, predicted values and quality criteria for the provided variables.
}
\details{
A typical predictor has the form dependent_Variable '~' independent_Variables.\cr The dependent_Variable has \strong{two} categories.\cr If there is more than one independent_Variable, they can be combined with a '+'.

The data frame must \strong{not} contain any missing values.\cr \strong{Metric} variables have to be of type \strong{numeric}. All \strong{other variables} have to be of type \strong{integer}.\cr The first variable in the dataset hat to be the dependent variable.\cr The scale of large numbers has to be reduced e.g. standardization.
}
\examples{
# generate Data
y <- as.integer(c(1,0,0,0,1,1,1,0,0,1))
x <- c(15,88,90,60,24,30,26,57,69,18)
z <- as.integer(c(3,2,2,1,3,3,2,1,1,3))
example_data <- data.frame(y,x,z)

# deploy QAS.func-Function
result <- QAS.func(y~x+z, data=example_data, weights=NULL, seed=NULL, tau = NULL)

}
\references{
King, G. & Zeng, L. (2001), Logistic Regression in Rare Events Data, \emph{Political Analysis}, No. 9 / 2001 \cr
Lipovetsky, S. (2014), Analytical closed-form solution for binary logit regression by categorical predictors, \emph{Journal of Applied Statistics}, No. 42 / 2015 \cr
Lipovetsky, S. & Conklin, M. (2014), Best-Worst Scaling in analytical closed-form solution, \emph{The Journal of Choice Modelling}, No. 10 / 2014 \cr
Stoltenberg, B. (2016), Using logit on big data - from iterative methods to analytical solutions, \emph{GfK Verein Working Paper Series}, No. 3 / 2016 \cr
}
\author{
This method is based on the research work of Stan Lipovetsky and Birgit Stoltenberg.
}
