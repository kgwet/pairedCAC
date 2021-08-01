#							paired t-test for agreement coefficients.r
#						 		  (February 20, 2021)
#Description: This script file contains a series of R functions for testing the difference between 2 agreement
#		 coefficients for statistical significance. The approach used here is based on the linearization method 
#    proposed Gwet (2016).A function is proposed for each agreement coefficient (Gwet's AC2, Fleiss' Kappa, 
#    Conger's Kappa, Krippendorff's Alpha, and Brennan-Prediger Coefficient.
#
#    Each function has only 2 required parameters, which are the 2 data frames containing the 2 sets of ratings for
#    the 2 groups of raters under comparison. Each input data file is in the form of nxr matrix or data frame showing 
#    the actual ratings each rater (column) assigned to each subject (in row). A typical table entry (i,g) represents 
#    the rating associated with subject i and rater g. 
#
#Author: Kilem L. Gwet, Ph.D.  (gwet@agreestat.com)
#
#
#     EXAMPLES OF SIMPLE CALLS OF THE 5 FUNCTIONS
#     ===========================================
# > ttest.fleiss(g1.ratings,g2.ratings)
# > ttest.ac2(g1.ratings,g2.ratings)
# > ttest.alpha(g1.ratings,g2.ratings)
# > ttest.conger(g1.ratings,g2.ratings)
# > ttest.bp(g1.ratings,g2.ratings)


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
#=====================================================================================
#ttest.ac2: This function performs the paired t-test for testing the difference  
#              between two correlated Gwet's AC2 coefficients for statistical significance. It
#              implements the linearization method of Gwet (2016).
#-------------
# This function takes 2 required parameters, which are the 2 groups of raters being compared:
# "g1.ratings" and "g2.ratings." Both datasets must have the exact same number of rows, and each
# column represents one rater and contains its ratings (numeric or alphabetic). 
# The user must exclude all subjects that are not rated by any rater.
# Bibliography:
# ------------
# Gwet, K. L. (2008). Computing inter-rater reliability and its variance in the presence of high agreement.
#                     British Journal of Mathematical and Statistical Psychology, 61, 29-48.
# Gwet, K. L. (2016). "Testing the Difference of Correlated Agreement Coefficients for Statistical Significance,
#                      Educational and Psychological Measurement, Vol 76(4) 609-637.
#======================================================================================
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
    invisible(c(ac2.coeff1, ac2.coeff2,coeff.diff,std.err,t.stat,p.value))
  }else{
    cat("Both datasets must have the same number of subjects. One has ",n1," subjects, while the other has ", n2," subjects.")
    invisible(-1)
  }
}

#=====================================================================================
#ttest.alpha: This function performs the paired t-test for testing the difference  
#              between two correlated Krippendorff's alpha coefficients for statistical significance. It
#              implements the linearization method of Gwet (2016).
#-------------
# This function takes 2 required parameters, which are the 2 groups of raters being compared:
# "g1.ratings" and "g2.ratings." Both datasets must have the exact same number of rows, and each
# column represents one rater and contains its ratings (numeric or alphabetic). 
# The user must exclude all subjects that are not rated by any rater.
# Bibliography:
# ------------
# Krippendorff, K. (1970). Estimating the reliability, systematic error, and random error of interval data.
#                          Educational and Psychological Measurement, 30, 61-70.
# Gwet, K. L. (2016). "Testing the Difference of Correlated Agreement Coefficients for Statistical Significance,
#                      Educational and Psychological Measurement, Vol 76(4) 609-637.
#======================================================================================
ttest.alpha <- function(g1.ratings,g2.ratings,weights="unweighted",conflev=0.95,N=Inf,print=TRUE){ 
  n2 <- nrow(g2.ratings)
  n1 <- nrow(g1.ratings)
  if (n2==n1){
    coeff2.i<-alpha.linear.i(g2.ratings,weights,conflev,N)
    coeff1.i<-alpha.linear.i(g1.ratings,weights,conflev,N)
    alpha.coeff1 = mean(coeff1.i)
    alpha.coeff2 = mean(coeff2.i)
    coeff.diff <- alpha.coeff2-alpha.coeff1
    n1i <- length(coeff1.i)
    n2i <- length(coeff2.i)
    if (n1i==n2i){
      di = coeff2.i-coeff1.i
      std.err <- sqrt(var(di)/n1)
      t.stat <- (alpha.coeff2-alpha.coeff1)/std.err
      p.value <- 2*(1-pt(abs(t.stat),n1-1))
    }else{
      std.err <- NA
      t.stat <- NA
      p.value <- NA
    }
    if (print==TRUE){
      cat("PAIRED T-TEST FOR TESTING THE DIFFERENCE BETWEEN 2 KRIPPENDORFF's Alpha COEFFICIENTS\n")
      cat("------------------------------------------------------------------------------------\n")
      cat("Alpha Coefficients: (Group 1: ", alpha.coeff1, ") -- (Group 2: ",alpha.coeff2,")\n")
      cat("Standard Error of the Differences: ", std.err,"\n")
      cat("Test Statistic: T= ",t.stat,"\n")
      cat("P-value: ",p.value)
    }
    invisible(c(alpha.coeff1, alpha.coeff2,coeff.diff,std.err,t.stat,p.value))
  }else{
    cat("Both datasets must have the same number of subjects. One has ",n1," subjects, while the other has ", n2," subjects.")
    invisible(-1)
  }
}

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

# ==============================================================
# This is an r function for trimming leading and trealing blanks
# ==============================================================
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x) }

#=====================================================================================
#fleiss.linear.i: This function produces all n linear components associated with  
#		           the n subjects, and the linearized version of the Fleiss' kappa coefficient.
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Fleiss, J. L. (1981). Statistical Methods for Rates and Proportions. John Wiley & Sons.
#======================================================================================
fleiss.linear.i <- function(ratings,weights="unweighted",conflev=0.95,N=Inf){ 
  ratings.mat <- as.matrix(ratings) 
  if (is.character(ratings.mat)){ratings.mat <- toupper(trim(ratings.mat))}
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 
  
  # creating a vector containing all categories used by the raters
  
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init)){
    categ <- sort(as.vector(na.omit(categ.init)))
  }else{
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- categ.init[nchar(categ.init)>0]
    categ <- sort(categ)
  }
  q <- length(categ)
  
  # creating the weights matrix
  
  if (is.character(weights)){
    if (weights=="quadratic")
      weights.mat<-quadratic.weights(categ)
    else if (weights=="ordinal")
      weights.mat<-ordinal.weights(categ)
    else if (weights=="linear")
      weights.mat<-linear.weights(categ)
    else if (weights=="radical")
      weights.mat<-radical.weights(categ)
    else if (weights=="ratio")
      weights.mat<-ratio.weights(categ)
    else if (weights=="circular")
      weights.mat<-circular.weights(categ)
    else if (weights=="bipolar")
      weights.mat<-bipolar.weights(categ)
    else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category
  
  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
    if (is.numeric(ratings.mat)){
      k.mis <-(ratings.mat==categ[k])
      in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
      agree.mat[,k] <- in.categ.k%*%rep(1,r) 
    }else{
      in.k <- (trim(ratings.mat)==categ[k])
      in.k[is.na(in.k)] <- FALSE 
      agree.mat[,k] <- in.k%*%rep(1,r)
    }
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
  
  # calculating fleiss's generalized kappa coefficient
  
  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more
  
  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat * (pi.vec%*%t(pi.vec)))
  fleiss.kappa <- (pa-pe)/(1-pe)
  
  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec
  
  pe.r2 <- pe*(ri.vec>=2)
  kappa.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  pi.vec.wk. <- weights.mat%*%pi.vec
  pi.vec.w.k <- t(weights.mat)%*%pi.vec
  pi.vec.w <- (pi.vec.wk. + pi.vec.w.k)/2
  
  pe.ivec <- (agree.mat%*%pi.vec.w)/ri.vec
  kappa.ivec.x <- kappa.ivec - 2*(1-fleiss.kappa) * (pe.ivec-pe)/(1-pe)
  invisible(kappa.ivec.x)
}
#===========================================================================================
#ac2.linear.i: This function produces all n linear components associated with  
#		           the n subjects, and the linearized version of the ac2 coefficient.
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Gwet, K. L. (2008). ``Computing inter-rater reliability and its variance in the presence of high
#		agreement." British Journal of Mathematical and Statistical Psychology, 61, 29-48.
#============================================================================================
ac2.linear.i <- function(ratings,weights="unweighted",conflev=0.95,N=Inf){ 
  ratings.mat <- as.matrix(ratings)
  if (is.character(ratings.mat)){ratings.mat <- toupper(trim(ratings.mat))}
  
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 
  
  # creating a vector containing all categories used by the raters
  
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init))
    categ <- sort(as.vector(na.omit(categ.init)))
  else{
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- categ.init[nchar(categ.init)>0]
    categ <- sort(categ)
  }
  q <- length(categ)
  
  # creating the weights matrix
  
  if (is.character(weights)){
    if (weights=="quadratic")
      weights.mat<-quadratic.weights(categ)
    else if (weights=="ordinal")
      weights.mat<-ordinal.weights(categ)
    else if (weights=="linear")
      weights.mat<-linear.weights(categ)
    else if (weights=="radical")
      weights.mat<-radical.weights(categ)
    else if (weights=="ratio")
      weights.mat<-ratio.weights(categ)
    else if (weights=="circular")
      weights.mat<-circular.weights(categ)
    else if (weights=="bipolar")
      weights.mat<-bipolar.weights(categ)
    else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category
  
  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
    if (is.numeric(ratings.mat)){
      k.mis <-(ratings.mat==categ[k])
      in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
      agree.mat[,k] <- in.categ.k%*%rep(1,r) 
    }else{
      in.k <- (trim(ratings.mat)==categ[k])
      in.k[is.na(in.k)] <- FALSE 
      agree.mat[,k] <- in.k%*%rep(1,r)
    }
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
  
  # calculating gwet's ac1 coefficient
  
  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more
  
  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  if (q>1){
    pe <- sum(weights.mat) * sum(pi.vec*(1-pi.vec)) / (q*(q-1))
  }else{pe <- 0}
  gwet.ac1 <- (pa-pe)/(1-pe)
  
  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec
  
  pe.r2 <- pe*(ri.vec>=2)
  ac1.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  if (q>1){
    pe.ivec <- (sum(weights.mat)/(q*(q-1))) * (agree.mat%*%(1-pi.vec))/ri.vec
  }else{
    pe.ivec <- matrix(0,n,1)
  }
  ac1.ivec.x <- ac1.ivec - 2*(1-gwet.ac1) * (pe.ivec-pe)/(1-pe)
  invisible(ac1.ivec.x)
}

#=====================================================================================
#alpha.linear.i: This function produces all n linear components associated with  
#		           the n subjects, and the linearized version of Krippendorff's alpha coefficient.
#-------------
#The algorithm used to compute krippendorff's alpha is very different from anything that was published on this topic. Instead,
#it follows the equations presented by K. Gwet (2014, p. 88)
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Gwet, K. (2014). Handbook of Inter-Rater Reliability: The Definitive Guide to Measuring the Extent of Agreement Among 
#	Multiple Raters, 4th Edition.  Advanced Analytics, LLC; 4th edition (September 10, 2014)
#Krippendorff (1970). "Bivariate agreement coefficients for reliability of data." Sociological Methodology,2,139-150
#Krippendorff (1980). Content analysis: An introduction to its methodology (2nd ed.), New-bury Park, CA: Sage.
#======================================================================================
alpha.linear.i <- function(ratings,weights="unweighted",conflev=0.95,N=Inf){ 
  ratings.mat <- as.matrix(ratings) 
  if (is.character(ratings.mat)){ratings.mat <- toupper(trim(ratings.mat))}
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 
  
  # creating a vector containing all categories used by the raters
  
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init)){
    categ <- sort(as.vector(na.omit(categ.init)))
  }else{
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- categ.init[nchar(categ.init)>0]
    categ <- sort(categ)
  }
  q <- length(categ)
  
  # creating the weights matrix
  
  if (is.character(weights)){
    if (weights=="quadratic")
      weights.mat<-quadratic.weights(categ)
    else if (weights=="ordinal")
      weights.mat<-ordinal.weights(categ)
    else if (weights=="linear")
      weights.mat<-linear.weights(categ)
    else if (weights=="radical")
      weights.mat<-radical.weights(categ)
    else if (weights=="ratio")
      weights.mat<-ratio.weights(categ)
    else if (weights=="circular")
      weights.mat<-circular.weights(categ)
    else if (weights=="bipolar")
      weights.mat<-bipolar.weights(categ)
    else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category
  
  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
    if (is.numeric(ratings.mat)){
      k.mis <-(ratings.mat==categ[k])
      in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
      agree.mat[,k] <- in.categ.k%*%rep(1,r) 
    }else{
      in.k <- (trim(ratings.mat)==categ[k])
      in.k[is.na(in.k)] <- FALSE 
      agree.mat[,k] <- in.k%*%rep(1,r)
    }
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
  
  # calculating krippendorff's alpha coefficient
  
  ri.vec <- agree.mat%*%rep(1,q)
  agree.mat<-as.matrix(agree.mat[(ri.vec>=2),])
  agree.mat.w <- agree.mat.w[(ri.vec>=2),]
  ri.vec <- ri.vec[(ri.vec>=2)]
  ri.mean <- mean(ri.vec)
  n <- nrow(agree.mat)
  epsi <- 1/sum(ri.vec)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%as.matrix(rep(1,q))
  pa <- (1-epsi)* sum(sum.q/(ri.mean*(ri.vec-1)))/n + epsi
  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/ri.mean))
  pe <- sum(weights.mat * (pi.vec%*%t(pi.vec)))
  krippen.alpha <- (pa-pe)/(1-pe)
  
  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  
  den.ivec <- ri.mean*(ri.vec-1)
  pa.ivec <- sum.q/den.ivec
  pa.v <- mean(pa.ivec)
  pa.ivec <- (1-epsi)*(pa.ivec-pa.v*(ri.vec-ri.mean)/ri.mean) + epsi
  krippen.ivec <- (pa.ivec-pe)/(1-pe)
  pi.vec.wk. <- weights.mat%*%pi.vec
  pi.vec.w.k <- t(weights.mat)%*%pi.vec
  
  pi.vec.w <- (pi.vec.wk. + pi.vec.w.k)/2
  
  pe.ivec <- (agree.mat%*%pi.vec.w)/ri.mean - sum(pi.vec) * (ri.vec-ri.mean)/ri.mean
  krippen.ivec.x <- krippen.ivec - 2*(1-krippen.alpha) * (pe.ivec-pe)/(1-pe)
  invisible(krippen.ivec.x)
}

#===========================================================================================
#conger.linear.i: This function produces all n linear components associated with  
#		           the n subjects, and the linearized version of Conger's kappa coefficient.
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Conger, A. J. (1980), ``Integration and Generalization of Kappas for Multiple Raters,"
#		Psychological Bulletin, 88, 322-328.
#======================================================================================
conger.linear.i <- function(ratings,weights="unweighted",conflev=0.95,N=Inf){ 
  ratings.mat <- as.matrix(ratings) 
  if (is.character(ratings.mat)){ratings.mat <- toupper(trim(ratings.mat))}
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 
  
  # creating a vector containing all categories used by the raters
  
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init))
    categ <- sort(as.vector(na.omit(categ.init)))
  else {
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- categ.init[nchar(categ.init)>0]
    categ <- sort(categ)
  }
  q <- length(categ)
  
  # creating the weights matrix
  
  if (is.character(weights)){
    if (weights=="quadratic")
      weights.mat<-quadratic.weights(categ)
    else if (weights=="ordinal")
      weights.mat<-ordinal.weights(categ)
    else if (weights=="linear")
      weights.mat<-linear.weights(categ)
    else if (weights=="radical")
      weights.mat<-radical.weights(categ)
    else if (weights=="ratio")
      weights.mat<-ratio.weights(categ)
    else if (weights=="circular")
      weights.mat<-circular.weights(categ)
    else if (weights=="bipolar")
      weights.mat<-bipolar.weights(categ)
    else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)

  # creating the nxq agreement matrix representing the distribution of raters by subjects and category
  
  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
    if (is.numeric(ratings.mat)){
      k.mis <-(ratings.mat==categ[k])
      in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
      agree.mat[,k] <- in.categ.k%*%rep(1,r) 
    }else{
      in.k <- (trim(ratings.mat)==categ[k])
      in.k[is.na(in.k)] <- FALSE 
      agree.mat[,k] <- in.k%*%rep(1,r)
    }
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
  
  # creating the rxq rater-category matrix representing the distribution of subjects by rater and category
  
  classif.mat <- matrix(0,nrow=r,ncol=q)
  for(k in 1:q){
    if (is.numeric(ratings.mat)){
      with.mis <-(t(ratings.mat)==categ[k])
      without.mis <- replace(with.mis,is.na(with.mis),FALSE)
      classif.mat[,k] <- without.mis%*%rep(1,n)
    }else{
      incat.k <- (t(ratings.mat)==categ[k])
      incat.k[is.na(incat.k)] <- FALSE
      classif.mat[,k] <- incat.k%*%rep(1,n)
    }
  }
  
  # calculating conger's kappa coefficient
  
  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more
  
  ng.vec <- classif.mat%*%rep(1,q)
  pgk.mat <- classif.mat/(ng.vec%*%rep(1,q))
  p.mean.k <- (t(pgk.mat)%*%rep(1,r))/r 
  s2kl.mat <- (t(pgk.mat)%*%pgk.mat - r * p.mean.k%*%t(p.mean.k))/(r-1)
  pe <- sum(weights.mat * (p.mean.k%*%t(p.mean.k) -  s2kl.mat/r))
  conger.kappa <- (pa-pe)/(1-pe)
  
  # calculating variance, stderr & p-value of conger's kappa coefficient
  
  bkl.mat <- (weights.mat+t(weights.mat))/2
  pe.ivec1 <- r*(agree.mat%*%t(t(p.mean.k)%*%bkl.mat))
  pe.ivec2 = rep(0,n)
  
  lamda.ig.mat=matrix(0,n,r)
  epsi.ig.mat <-1-is.na(ratings.mat)
  epsi.ig.mat <- replace(epsi.ig.mat,is.na(epsi.ig.mat),FALSE)
  for(k in 1:q){
    lamda.ig.kmat=matrix(0,n,r)
    for(l in 1:q){
      delta.ig.mat <- (ratings.mat==categ[l])
      delta.ig.mat <- replace(delta.ig.mat,is.na(delta.ig.mat),FALSE)
      lamda.ig.kmat <- lamda.ig.kmat + weights.mat[k,l] * (delta.ig.mat - (epsi.ig.mat - rep(1,n)%*%t(ng.vec/n)) * (rep(1,n)%*%t(pgk.mat[,l])))
    }
    lamda.ig.kmat = lamda.ig.kmat*(rep(1,n)%*%t(n/ng.vec))
    lamda.ig.mat = lamda.ig.mat+ lamda.ig.kmat*(r*mean(pgk.mat[,k]) - rep(1,n)%*%t(pgk.mat[,k]))
  }
  pe.ivec <- (lamda.ig.mat%*%rep(1,r)) / (r*(r-1))
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec
  pe.r2 <- pe*(ri.vec>=2)
  conger.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe) 
  conger.ivec.x <- conger.ivec - 2*(1-conger.kappa) * (pe.ivec-pe)/(1-pe)
  invisible(conger.ivec.x)
}
#===========================================================================================
#bp.linear.i: This function produces all n linear components associated with  
#		           the n subjects, and the linearized version of B-P's coefficient.
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Brennan, R.L., and Prediger, D. J. (1981). ``Coefficient Kappa: some uses, misuses, and alternatives."
#           Educational and Psychological Measurement, 41, 687-699.
#======================================================================================
bp.linear.i <- function(ratings,weights="unweighted",conflev=0.95,N=Inf){ 
  ratings.mat <- as.matrix(ratings) 
  if (is.character(ratings.mat)){ratings.mat <- toupper(trim(ratings.mat))}
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 
  
  # creating a vector containing all categories used by the raters
  
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init))
    categ <- sort(as.vector(na.omit(categ.init)))
  else{
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- categ.init[nchar(categ.init)>0]
    categ <- sort(categ)
  }
  q <- length(categ)
  
  # creating the weights matrix
  
  if (is.character(weights)){
    if (weights=="quadratic")
      weights.mat<-quadratic.weights(categ)
    else if (weights=="ordinal")
      weights.mat<-ordinal.weights(categ)
    else if (weights=="linear")
      weights.mat<-linear.weights(categ)
    else if (weights=="radical")
      weights.mat<-radical.weights(categ)
    else if (weights=="ratio")
      weights.mat<-ratio.weights(categ)
    else if (weights=="circular")
      weights.mat<-circular.weights(categ)
    else if (weights=="bipolar")
      weights.mat<-bipolar.weights(categ)
    else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category
  
  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
    if (is.numeric(ratings.mat)){
      k.mis <-(ratings.mat==categ[k])
      in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
      agree.mat[,k] <- in.categ.k%*%rep(1,r) 
    }else{
      in.k <- (trim(ratings.mat)==categ[k])
      in.k[is.na(in.k)] <- FALSE 
      agree.mat[,k] <- in.k%*%rep(1,r)
    }
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
  
  # calculating gwet's ac1 coefficient
  
  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more
  
  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  if (q>1){
    pe <- sum(weights.mat) / (q^2)
  }else{pe <- 0}
  bp.coeff <- (pa-pe)/(1-pe)
  
  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec
  
  pe.r2 <- pe*(ri.vec>=2)
  bp.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  invisible(bp.ivec)
}
