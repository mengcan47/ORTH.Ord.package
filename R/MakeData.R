
# Module: MakeData
#
# to create analytical dataset
#
# @param sim_y is column name for response variable y
# @param sim_X is a vector of column names for covariates X
# @param sim_Z is a vector of the labels for alpha
# @param sim_data is the dataset with y and X
# @param sim_zdata is the dataset for pairewise correlation
# @param sim_cluster is the column name for cluster indicator
# @param ordinal_L is the number of levels of ordinal response minus 1
# @param independence is an indicator for whether independence working correlation is used
# @param sim_zcutpoint is an indicator for whether correlation parameters are dependent on cutpoint
# @param sim_xpartial is an indicator


MakeData<-function(sim_y,sim_X,sim_Z,sim_data,sim_zdata,sim_cluster,ordinal_L,independence,sim_zcutpoint=0,sim_xpartial=NULL){
  mb<-sim_data[order(sim_data[[sim_cluster]]),]
  y<-c(mb[,paste0(sim_y)])
  if (length(sim_X)!=0){X<-data.frame(mb[,paste0(sim_X)])}
  Sample_N<-data.frame(table(mb[,sim_cluster]))##_sample dataset in SAS
  colnames(Sample_N)<-c("cluster","N")
  n<-c(Sample_N$N)#size of each clusters
  L<-ordinal_L
  id<-c(mb[,sim_cluster])

  if (length(sim_xpartial)!=0) {
    Xp<-data.frame(mb[,paste0(sim_xpartial)])
  }
  ###Expanding y, X, and id for ordinal level indicators;
  yc<-y
  if (length(sim_X)!=0){xc<-X}

  idc<-id

  X<-NULL
  y<-NULL
  id<-NULL

  y<-ExpandY(yc)

  for (i in 1:length(idc)){
    idi<-kronecker(rep(1,L),idc[i])
    if (length(sim_X)!=0){
      Xi<-cbind(diag(L),kronecker(rep(1,L),data.matrix(xc[i,])))
    } else {
      Xi<-diag(L)
    }
    ###Non-porportional covariates;
    if (length(sim_xpartial)!=0) {
      Xi<-cbind(Xi,kronecker(data.matrix(Xp[i,]),c(rep(0,(L-1)),diag(L-1))))
    }
    id<-c(id,idi)
    X<-rbind(X,Xi)
  }


  ###Z matrix;
  Z<-NULL
  if (independence==0){
    Zc<-sim_zdata[,c(paste0(sim_cluster),paste0(sim_Z))]
    for (i in 1:length(n)) {
      c_i<-rep(c(0:(L-1)),n[i])
      z_i_temp<-data.matrix(Zc[which(Zc[,1]==i),])
      colnames(z_i_temp)<-NULL
      rownames(z_i_temp)<-NULL
      z_i<-kronecker(z_i_temp,rep(1,L^2))
      z_i_temp<-NULL
      if (sim_zcutpoint==1){
        A<-rbind(t(rep(0,(L^2-1))),diag(L^2-1))
        z_i_temp<-data.matrix(Zc[which(Zc[,1]==i),2:ncol(Zc)])
        z_i<-cbind(z_i,kronecker(z_i_temp,A))
        z_i_temp<-NULL
      }
      Z<-rbind(Z,z_i)
    }
  }

  Result_data<-list(y,X,Z,id,n,L)
  return(Result_data)
}
