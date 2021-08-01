circular.weights<-function(categ){
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
  U = xmax-xmin+1
  for(k in 1:q){
    for(l in 1:q){
      weights[k,l] <- (sin(pi*(categ.vec[k]-categ.vec[l])/U))^2
    }
  }
  weights <- 1-weights/max(weights)
  return (weights)
}
