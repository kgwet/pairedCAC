#' @title Paired t-test for the difference of 2 Brennan-Prediger coefficients
#'
#' @description
#' This function performs the paired t-test for testing the difference between
#' two correlated Brennan-Prediger coefficients for statistical significance.
#' It implements the linearization method of Gwet (2016).
#'
#' @param g1.ratings is a mandatory parameter representing the first data frame
#' of ratings.
#' @param g2.ratings is a mandatory parameter representing the second data frame
#' of ratings.
#' @param weights is an optional parameter that defines the weights needed in a
#' weighted analysis. It's default value is ``unweighted'' which requests the
#' unweighted analysis.
#' @param conflev, is an optional parameter representing the confidence level.
#' It's default value is 0.95.
#' @param N=, is an optional parameter representing the size of the subject
#' population. It's default value is infinity.
#'
#' @details
#' This function takes 2 required parameters, which are the 2 groups of raters
#' being compared: "g1.ratings" and "g2.ratings." Both datasets must have the
#' exact same number of rows, and each column represents one rater and
#' contains its ratings (numeric or alphabetic).
#' The user must exclude all subjects that are not rated by any rater.
#'
#' @references
#' Brennan, R.L. & Prediger, D.J. (1981). ``Coefficient kappa: some uses,
#' misuses, and alternatives. \emph{Educational and Psychological Measurement},
#' 41, 687-699.
#'
#' Gwet, K.L. (2016). ``Testing the Difference of Correlated Agreement
#' Coefficients for Statistical Significance,'' \emph{Educational and
#' Psychological Measurement}, \strong{76}(4) 609-637.
#' @export
ttest.bp <- function(g1.ratings,g2.ratings,weights="unweighted",conflev=0.95,N=Inf){
  n2 <- nrow(g2.ratings)
  n1 <- nrow(g1.ratings)
  if (n2==n1){
    coeff2.i<-bp.linear.i(g2.ratings,weights,conflev,N)
    coeff1.i<-bp.linear.i(g1.ratings,weights,conflev,N)
    di = coeff2.i$icoeff-coeff1.i$icoeff
    bp.coeff1 = mean(coeff1.i$icoeff)
    bp.coeff2 = mean(coeff2.i$icoeff)
    coeff.diff <- bp.coeff2-bp.coeff1
    std.err <- sqrt(var(di)/n1)
    t.stat <- (bp.coeff2-bp.coeff1)/std.err
    p.value <- 2*(1-pt(abs(t.stat),n1-1))
    n.obs <- n2
    n.raters1 <- ncol(g1.ratings)
    n.raters2 <- ncol(g2.ratings)
    weight.mat <- coeff1.i$weights
    df.out <- data.frame(bp.coeff1, bp.coeff2,coeff.diff,std.err,t.stat,p.value,n.obs,n.raters1,n.raters2)
  }else{
    cat("Both datasets must have the same number of subjects. One has ",n1," subjects, while the other has ", n2," subjects.")
    df.out <- NULL
    weight.mat <- NULL
  }
  return(list("test"=df.out,"weights"=weight.mat))
}
