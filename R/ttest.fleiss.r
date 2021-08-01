#=====================================================================================
#ttest.fleiss: This function performs the paired t-test for testing the difference
#              between two correlated Fleiss kappa coefficients for statistical significance. It
#              implements the linearization method of Gwet (2016).
#-------------
# This function takes 2 required parameters, which are the 2 groups of raters being compared:
# "g1.ratings" and "g2.ratings." Both datasets must have the exact same number of rows, and each
# column represents one rater and contains its ratings (numeric or alphabetic).
# The user must exclude all subjects that are not rated by any rater.
# Bibliography:
# ------------
# Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters.
#                        Psychological Bulletin, 76, 378-382.
# Gwet, K. L. (2016). "Testing the Difference of Correlated Agreement Coefficients for Statistical Significance,
#                      Educational and Psychological Measurement, Vol 76(4) 609-637.
#======================================================================================
ttest.fleiss <- function(g1.ratings,g2.ratings,weights.p="unweighted",conflev=0.95,N=Inf,print=TRUE){
  n2 <- nrow(g2.ratings)
  n1 <- nrow(g1.ratings)
  if (n2==n1){
    coeff2.i<-fleiss.linear.i(g2.ratings,weights=weights.p,conflev,N)
    coeff1.i<-fleiss.linear.i(g1.ratings,weights=weights.p,conflev,N)
    di = coeff2.i-coeff1.i
    fleiss.coeff1 = mean(coeff1.i)
    fleiss.coeff2 = mean(coeff2.i)
    coeff.diff <- fleiss.coeff2-fleiss.coeff1
    std.err <- sqrt(var(di)/n1)
    t.stat <- (fleiss.coeff2-fleiss.coeff1)/std.err
    p.value <- 2*(1-pt(abs(t.stat),n1-1))
    if (print==TRUE){
      cat("PAIRED T-TEST FOR TESTING THE DIFFERENCE BETWEEN 2 FLEISS AGREEMENT COEFFICIENTS\n")
      cat("--------------------------------------------------------------------------------\n")
      cat("Fleiss Kappa Coefficients: (Group 1: ", fleiss.coeff1, ") -- (Group 2: ",fleiss.coeff2,")\n")
      cat("Standard Error of the Differences: ", std.err,"\n")
      cat("Test Statistic: T= ",t.stat,"\n")
      cat("P-value: ",p.value)
    }
    invisible(c(fleiss.coeff1, fleiss.coeff2,coeff.diff,std.err,t.stat,p.value))
  }else{
    cat("Both datasets must have the same number of subjects. One has ",n1," subjects, while the other has ", n2," subjects.")
    invisible(-1)
  }
}
