ordinal.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  categ.vec<-1:length(categ)
  for(k in 1:q){
    for(l in 1:q){
      nkl <- max(k,l)-min(k,l)+1
      weights[k,l] <- nkl * (nkl-1)/2
    }
  }
  weights <- 1-weights/max(weights)
  return (weights)
}
