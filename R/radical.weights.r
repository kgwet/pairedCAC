radical.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  if (is.numeric(categ)) {
    categ.vec <- sort(categ)
  }
  else {
    categ.vec<-1:length(categ)
  }
  xmin<-min(categ.vec)
  xmax<-max(categ.vec)
  for(k in 1:q){
    for(l in 1:q){
      weights[k,l] <- 1-sqrt(abs(categ.vec[k]-categ.vec[l]))/sqrt(abs(xmax-xmin))
    }
  }
  return (weights)
}
