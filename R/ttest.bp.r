#=====================================================================================
#ttest.bp: This function performs the paired t-test for testing the difference
#              between two correlated Brennan-Prediger coefficients for statistical significance. It
#              implements the linearization method of Gwet (2016).
#-------------
# This function takes 2 required parameters, which are the 2 groups of raters being compared:
# "g1.ratings" and "g2.ratings." Both datasets must have the exact same number of rows, and each
# column represents one rater and contains its ratings (numeric or alphabetic).
# The user must exclude all subjects that are not rated by any rater.
# Bibliography:
# ------------
# Brennan, R.L. & Prediger, D. J. (1981). Coefficient kappa: some uses, misuses, and alternatives.
#                       Educational and Psychological Measurement, 41, 687-699.
# Gwet, K. L. (2016). "Testing the Difference of Correlated Agreement Coefficients for Statistical Significance,
#                      Educational and Psychological Measurement, Vol 76(4) 609-637.
#======================================================================================
ttest.bp <- function(g1.ratings,g2.ratings,weights="unweighted",conflev=0.95,N=Inf,print=TRUE){
  n2 <- nrow(g2.ratings)
  n1 <- nrow(g1.ratings)
  if (n2==n1){
    coeff2.i<-bp.linear.i(g2.ratings,weights,conflev,N)
    coeff1.i<-bp.linear.i(g1.ratings,weights,conflev,N)
    di = coeff2.i-coeff1.i
    bp.coeff1 = mean(coeff1.i)
    bp.coeff2 = mean(coeff2.i)
    coeff.diff <- bp.coeff2-bp.coeff1
    std.err <- sqrt(var(di)/n1)
    t.stat <- (bp.coeff2-bp.coeff1)/std.err
    p.value <- 2*(1-pt(abs(t.stat),n1-1))
    if (print==TRUE){
      cat("PAIRED T-TEST FOR TESTING THE DIFFERENCE BETWEEN 2 BRENNAN-PREDIGER COEFFICIENTS\n")
      cat("--------------------------------------------------------------------------------\n")
      cat("Brennan-Prediger Coefficients: (Group 1: ", bp.coeff1, ") -- (Group 2: ",bp.coeff2,")\n")
      cat("Standard Error of the Differences: ", std.err,"\n")
      cat("Test Statistic: T= ",t.stat,"\n")
      cat("P-value: ",p.value)
    }
    invisible(c(bp.coeff1, bp.coeff2,coeff.diff,std.err,t.stat,p.value))
  }else{
    cat("Both datasets must have the same number of subjects. One has ",n1," subjects, while the other has ", n2," subjects.")
    invisible(-1)
  }
}
