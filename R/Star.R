
# Module: Star
#
# Modules to manage conditional residuals for ordinal variables
#
#
# @param  y is a vector generated from ExpandY function
# @param  L is \code{max}(y), y=0,1,...,L.
#
# @return a vector Y*, see Heagerty and Zeger (1996).
#


Star<-function(y,L){
  obs<-length(y)/L
  for (j in 0:(obs-2)) {
    a<-kronecker(data.matrix(y[(j*L+1):((j+1)*L)]),rep(1,L))
    b<-rep(a,(obs-j-1))
    if (j==0){o<-b}else{o<-c(o,b)}
  }
  return(o)
}
