
# Module: FirstBeta
#
# to initialize beta parameter
#
# @param ibeta is default to NULL if not specified
# @param y is response vector
# @param X is design matrix
#
# @return initial value of beta

firstbeta<-function(ibeta=NULL,y,X){
  if (length(ibeta)!=0){
    beta<-ibeta
  }else{
    z<-y+y-1
    beta<-solve(t(X)%*%X,t(X)%*%z)
    for (i in 1:2) {
      u<-1/(1+exp(-X%*%beta))
      v<-u*(1-u)
      z<-t(X)%*%(y-u)
      d<-solve(t(X)%*%(X*c(v)),z)
      beta<-beta+d
    }
  }
  return(beta)
}
