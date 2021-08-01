#=====================================================================================
#ttest.conger: This function performs the paired t-test for testing the difference
#              between two correlated Conger's Kappa coefficients for statistical significance. It
#              implements the linearization method of Gwet (2016).
#-------------
# This function takes 2 required parameters, which are the 2 groups of raters being compared:
# "g1.ratings" and "g2.ratings." Both datasets must have the exact same number of rows, and each
# column represents one rater and contains its ratings (numeric or alphabetic).
# The user must exclude all subjects that are not rated by any rater.
# Bibliography:
# ------------
# Conger, A. J. (1980). Integration and generalization of kappas for multiple raters.
#                       Psychological Bulletin, 88, 322-328.
# Gwet, K. L. (2016). "Testing the Difference of Correlated Agreement Coefficients for Statistical Significance,
#                      Educational and Psychological Measurement, Vol 76(4) 609-637.
#======================================================================================
ttest.conger <- function(g1.ratings,g2.ratings,weights="unweighted",conflev=0.95,N=Inf,print=TRUE){
  n2 <- nrow(g2.ratings)
  n1 <- nrow(g1.ratings)
  if (n2==n1){
    coeff1.i<-conger.linear.i(g1.ratings,weights,conflev,N)
    coeff2.i<-conger.linear.i(g2.ratings,weights,conflev,N)
    di = coeff2.i-coeff1.i
    conger.coeff1 = mean(coeff1.i)
    conger.coeff2 = mean(coeff2.i)

    coeff.diff <- conger.coeff2-conger.coeff1
    std.err <- sqrt(var(di)/n1)
    t.stat <- (conger.coeff2-conger.coeff1)/std.err
    p.value <- 2*(1-pt(abs(t.stat),n1-1))
    if (print==TRUE){
      cat("PAIRED T-TEST FOR TESTING THE DIFFERENCE BETWEEN 2 CONGER's KAPPA COEFFICIENTS\n")
      cat("------------------------------------------------------------------------------\n")
      cat("COnger's Coefficients: (Group 1: ", conger.coeff1, ") -- (Group 2: ",conger.coeff2,")\n")
      cat("Standard Error of the Differences: ", std.err,"\n")
      cat("Test Statistic: T= ",t.stat,"\n")
      cat("P-value: ",p.value)
    }
    invisible(c(conger.coeff1, conger.coeff2,coeff.diff,std.err,t.stat,p.value))
  }else{
    cat("Both datasets must have the same number of subjects. One has ",n1," subjects, while the other has ", n2," subjects.")
    invisible(-1)
  }
}
