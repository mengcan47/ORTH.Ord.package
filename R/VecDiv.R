
# Module: VecDiv
# @description This function will perform vector division by column
# @param M: A n by m matrix
# @param v: A n by 1 vector
#
# @return A matrix whose column i equals M[,i] divided by v
# @examples
# \dontrun{
# a<-matrix(c(1,2,3,4,5,6),3,2,byrow = T)
# b<-c(2,2,2)
# VecDiv(a,b)}


VecDiv<-function(M,v){
  n_c<-ncol(M)
  A<-matrix(0,nrow(M),n_c)
  for (i in 1:n_c){
    A[,i]<-M[,i]/v;
  }
  return(A)
}


