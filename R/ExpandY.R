
# Module:  ExpandY
#
# To format ordinal response Y as a vector of binary indicators for each observation
#
# @return a vector consist of binary indicators
# @param y is a vector


ExpandY<-function(y){
  L<-max(y)
  e<-ifelse(y==0,1,0)
  for (i in 1:(L-1)) {
    e<-cbind(e,ifelse(y==i,1,0))
  }
  d<-c(t(e))
  return(d)
}
