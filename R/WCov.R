
# Module: WCov
#
# Creates covariance matrix for beta estimating equation for clustered ordinal outcomes.
#
# @param A is variance of Y_i within each observation
# @param mujk is a vector of expected value of the product Y_ij and Y_ik
# @param mu is a vector of expected value of Y_ij
# @param L is \code{max}(y), y=0,1,...,L.
# @return covariance matrix for beta estimating equation


WCov<-function(A,mujk,mu,L,invT){
  #Note the ordering of the vector mujk is not the same as the order of elements in the Cov(Y_i) matrix.;
  ni<-nrow(mu)/L
  V<-A
  row_n<-L+1
  col_n<-1
  j<-1

  for (k in 1:(ni*(ni-1)/2*L^2)){
    V[row_n,col_n]<-mujk[k]-mu[row_n]*mu[col_n]
    V[col_n,row_n]<-V[row_n,col_n]
    if ((k%%L)==0){
      col_n<-col_n+1
      row_n<-row_n-L}
    if ((k%%(L^2))==0){
      row_n<-row_n+L
      col_n<-col_n-L
    }
    if (row_n==ni*L){
      j<-j+1
      row_n<-j*L
      col_n<-col_n+L
    }
    row_n<-row_n+1
  }
  #Note also that to estimate alpha, cumulative indicators were used, so that Cov(Y_i) has to be adjusted.;
  for (i in 1:(ni-1)){
    for (j in i:(ni-1)) {
      B<-V[(j*L+1):((j+1)*L),((i-1)*L+1):(i*L)]
      V[(j*L+1):((j+1)*L),((i-1)*L+1):(i*L)]<-invT%*%B%*%t(invT)
      V[((i-1)*L+1):(i*L),(j*L+1):((j+1)*L)]<-t(invT%*%B%*%t(invT))
    }
  }

  return(V)
}
