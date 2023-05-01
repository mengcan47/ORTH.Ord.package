
# Module: RCov
#
# to calculate covariance matrix P_i for correlation estimating equation
#
# @return element of covariance matrix for alpha estimating equation


Rcov<-function(mujk,muj,muk,b1,b2,L,rrow,ccol){
  i_star<-kronecker(c(1:L),rep(1,L))
  i_dstar<-kronecker(rep(1,L),c(1:L))
  i_jk<-cbind(i_star,i_dstar)

  c1<-i_star[rrow]
  c2<-i_dstar[rrow]
  c3<-i_star[ccol]
  c4<-i_dstar[ccol]
  min_c1c3<-min(c1,c3)
  min_c2c4<-min(c2,c4)

  tb<-which(i_jk[,1]==min_c1c3&i_jk[,2]==min_c2c4)

  ###Preissers;

  A<-(mujk[tb]-b1[rrow]*mujk[which(i_jk[,1]==min_c1c3&i_jk[,2]==c4)]
      -b1[ccol]*mujk[which(i_jk[,1]==min_c1c3&i_jk[,2]==c2)]
      -b2[rrow]*mujk[which(i_jk[,1]==c3&i_jk[,2]==min_c2c4)]
      -b2[ccol]*mujk[which(i_jk[,1]==c1&i_jk[,2]==min_c2c4)]
      +(muj[ccol]*b1[ccol] + muk[ccol]*b2[ccol]-mujk[ccol])*mujk[rrow]
      +(muj[rrow]*b1[rrow] + muk[rrow]*b2[rrow]-mujk[rrow])*mujk[ccol]
      +b1[rrow]*b2[ccol]*mujk[which(i_jk[,1]==c1&i_jk[,2]==c4)]
      +b1[ccol]*b2[rrow]*mujk[which(i_jk[,1]==c3&i_jk[,2]==c2)]
      +b1[rrow]*b1[ccol]*muj[tb]+b2[rrow]*b2[ccol]*muk[tb]
      -(muj[rrow]*b1[rrow] + muk[rrow]*b2[rrow])*(muj[ccol]*b1[ccol] + muk[ccol]*b2[ccol])
      +mujk[rrow]*mujk[ccol]
  )
  return(A)
}
