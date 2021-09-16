#' @title Create a weight matrix based raw data and weight set name
#'
#' @description The weight.mat.fn function produces the weight matrix needed
#' to compute weighted agreement coefficients for categorical ratings.
#' @param ratings is a mandatory parameter representing the input data frame
#' of raw ratings.
#' @param weights is an optional parameter that defines the weights needed in a
#' weighted analysis. It's default value is ``unweighted'' which requests the
#' unweighted analysis.
#'
#' @references
#' Gwet, K. L. (2021). \emph{Handbook of Inter-Rater Reliability - Volume 1: Analysis
#' of Categorical Ratings}. AgreeStat Analytics.
#' @export
weight.mat.fn <- function(ratings,weights="unweighted"){
  # creating the weights matrix based on the weight description
  ratings.mat <- as.matrix(ratings)
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init)){
    categ <- sort(as.vector(stats::na.omit(categ.init)))
  }else{
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- categ.init[nchar(categ.init)>0]
    categ <- sort(categ)
  }
  q <- length(categ)
  if (is.character(weights)){
    w.name <- weights
    if (weights=="quadratic") weights.mat<-quadratic.weights(categ)
    else if (weights=="ordinal") weights.mat<-ordinal.weights(categ)
    else if (weights=="linear") weights.mat<-linear.weights(categ)
    else if (weights=="radical") weights.mat<-radical.weights(categ)
    else if (weights=="ratio") weights.mat<-ratio.weights(categ)
    else if (weights=="circular") weights.mat<-circular.weights(categ)
    else if (weights=="bipolar") weights.mat<-bipolar.weights(categ)
    else weights.mat<-identity.weights(categ)
  }else{
    w.name <- "Custom weights"
    weights.mat= as.matrix(weights)
    if (sum(weights.mat)==q) w.name <- "Unweighted"
  }
  return(weights.mat)
}
