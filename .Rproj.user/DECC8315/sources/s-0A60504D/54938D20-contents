#' @title Paired t-test for the difference of 2 AC2 coefficients.
#'
#' @description The ttest.ac2 function performs the paired t-test for testing
#' the difference between two correlated Gwet's \eqn{AC_2} coefficients for
#' statistical significance. It implements the linearization method of
#' Gwet (2016).
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
#' @details The 2 required parameters "g1.ratings" and "g2.ratings",
#' which are the 2 groups of raters being compared, must have the
#' exact same number of rows, and each column represents one rater and
#' contains its ratings (numeric or alphabetic). All subjects that are not
#' rated by any rater must be excluded from the dataset.
#'
#' @references
#' Gwet, K. L. (2008). Computing inter-rater reliability and its variance in
#' the presence of high agreement. \emph{British Journal of Mathematical and
#' Statistical Psychology}, 61, 29-48.
#'
#' Gwet, K. L. (2016). Testing the Difference of Correlated Agreement
#' Coefficients for Statistical Significance, \emph{Educational and
#' Psychological Measurement}, Vol 76(4) 609-637.
ttest.ac2 <- function(g1.ratings,g2.ratings,weights="unweighted",conflev=0.95,N=Inf,print=TRUE){
  n2 <- nrow(g2.ratings)
  n1 <- nrow(g1.ratings)
  if (n2==n1){
    coeff2.i<-ac2.linear.i(g2.ratings,weights,conflev,N)
    coeff1.i<-ac2.linear.i(g1.ratings,weights,conflev,N)
    di = coeff2.i-coeff1.i
    ac2.coeff1 = mean(coeff1.i)
    ac2.coeff2 = mean(coeff2.i)
    coeff.diff <- ac2.coeff2-ac2.coeff1
    std.err <- sqrt(var(di)/n1)
    t.stat <- (ac2.coeff2-ac2.coeff1)/std.err
    p.value <- 2*(1-pt(abs(t.stat),n1-1))
    if (print==TRUE){
      cat("PAIRED T-TEST FOR TESTING THE DIFFERENCE BETWEEN 2 GWET's AC2 AGREEMENT COEFFICIENTS\n")
      cat("------------------------------------------------------------------------------------\n")
      cat("AC2 Coefficients: (Group 1: ", ac2.coeff1, ") -- (Group 2: ",ac2.coeff2,")\n")
      cat("Standard Error of the Differences: ", std.err,"\n")
      cat("Test Statistic: T= ",t.stat,"\n")
      cat("P-value: ",p.value)
    }
    df.out <- data.frame(ac2.coeff1, ac2.coeff2,coeff.diff,std.err,t.stat,p.value)
    #invisible(c(ac2.coeff1, ac2.coeff2,coeff.diff,std.err,t.stat,p.value))
  }else{
    cat("Both datasets must have the same number of subjects. One has ",n1," subjects, while the other has ", n2," subjects.")
    #invisible(-1)
    df.out <- NULL
  }
  return(df.out)
}
