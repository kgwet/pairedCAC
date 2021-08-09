#===========================================================================================
#pa.linear.i: This function produces all n linear components of the percent
#   agreement coefficient, associated with the n subjects.
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. You should exclude all subjects that are not
#rated by any rater.
#Bibliography:
#Gwet, K.L. (2016). Handbook of Inter-Rater Reliability, 5th Edition.
#     Volume 1: Analysis of Categorical Ratings. AgreeStat Analytics.
#======================================================================================
pa.linear.i <- function(ratings,weights="unweighted",conflev=0.95,N=Inf){
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

  # calculating the n linear components of the percent agreement pa

  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)

  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec
  pa.ivec <- (n/n2more)*pa.ivec
  return(list("icoeff"=pa.ivec,"weights"=weights.mat))
}
