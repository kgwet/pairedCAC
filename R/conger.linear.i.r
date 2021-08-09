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
  return(list("icoeff"=conger.ivec.x,"weights"=weights.mat))
}
