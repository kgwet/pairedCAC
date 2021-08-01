bipolar.weights<-function(categ){
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
      if (k!=l)
        weights[k,l] <- (categ.vec[k]-categ.vec[l])^2 / (((categ.vec[k]+categ.vec[l])-2*xmin)*(2*xmax-(categ.vec[k]+categ.vec[l])))
      else weights[k,l] <- 0
    }
  }
  weights <- 1-weights/max(weights)
  return (weights)
}
