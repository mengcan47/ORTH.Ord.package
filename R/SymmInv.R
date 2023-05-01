
#########################################################################################
#
# MODULE: SymmInv
# Returns an inverse of a symmetric matrix
#
# INPUT
# A: Symmetric matrix
#
# OUTPUT
# inv(A)
#
#########################################################################################

SymmInv<- function(A) {
  A<-as.matrix(A)
  if( det(A) > 0.0001) {
    B<-chol(A)
    C<-MASS::ginv(B)
    ainv<-C %*% t(C)
  }
  else {
    ainv = MASS::ginv(A)
  }
  return(ainv)
}
