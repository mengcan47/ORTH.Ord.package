
#######################################################################
#Module: Score
#
#Input: beta, alpha
#Output:
#U- Score vector
#UUt- Sum of U_i*U_i` across all clusters
#DVD- Approximate information matrix
#######################################################################


Score<-function(beta,alpha,y,X,Z,L,n,id,p,q,sim_Z,omega,MMORTH,independence){

  if (independence==0){
    U<-rep(0,p+q)
    DVD<-matrix(0L,nrow=p+q,ncol=p+q)
    UUt<-DVD
  } else {
    U<-rep(0,p)
    DVD<-matrix(0L,nrow=p,ncol=p)
    UUt<-DVD
  }


  ###Cumulative transform matrix ;
  T<-matrix(1L,nrow=L,ncol=L)
  for (row_n in 1:L) {
    for (col_n in 1:L) {
      if (col_n>row_n){T[row_n,col_n]<-0}
    }
  }
  invT<-solve(T)

  #record elements for bias correction: D, invV, Da, iVa
  D_list<-NULL
  iV_list<-NULL
  Da_list<-NULL
  iVa_list<-NULL
  mu_list<-NULL
  R_list<-NULL

  for (i in 1:length(n)) {
    D<-NULL
    obs<-NULL
    nobs<-NULL
    y_i<-NULL
    X_i<-NULL
    A_i<-NULL
    m_i<-NULL
    p_i<-NULL
    Tm_i<-NULL#not appear in the code below


    obs<-which(id==i)
    nobs<-length(obs)
    y_i<-y[obs]
    X_i<-X[obs,]

    if (independence==0){
      U_i<-rep(0,p+q)
      DVD_i<-matrix(0L,nrow = p+q,ncol = p+q)
    } else {
      U_i<-rep(0,p)
      DVD_i<-matrix(0L,nrow = p,ncol = p)
    }


    for (j in 1:n[i]) {
      X_it<-X_i[((j-1)*L+1):(j*L),]
      m_it<-1/(1+exp(-X_it%*%beta))
      p_it<-invT%*%m_it
      D_t<-t(t(X_it)%*%diag(c(T%*%p_it*(1-T%*%p_it)))%*%t(invT))
      m_i<-rbind(m_i,m_it)
      p_i<-rbind(p_i,p_it)

      if (j==1){D<-D_t}else{D<-rbind(D,D_t)}

      A_it<-diag(c(p_it))-p_it%*%t(p_it)
      #A_it<-diag(c(p_it*(1-p_it))) #--modified 20220121

      if (j==1){A_i<-A_it}else{A_i<-magic::adiag(A_i,A_it)}
    }


    if (independence==1 | n[i]==1){
      V_i<-A_i
    }else{###For alpha and working covariance estimation;
      mY_i<-NULL
      z_i<-NULL
      for (j in 1:n[i]) {###Cummulative indicators used to model odds ratio;
        mY_it<-T%*%y_i[((j-1)*L+1):(j*L)]
        mY_i<-rbind(mY_i,mY_it)
      }

      Ystar<-Star(mY_i,L)
      Ydstar<-DStar(mY_i,L)

      z_i<-Z[which(Z[,1]==i),2:ncol(Z)]

      gamma<-z_i%*%alpha #Log pairwise OR
      muj<-Star(m_i,L)
      muk<-DStar(m_i,L)
      mujk<-VecMardia(muj,muk,gamma)

      ###Working covariance for Beta estimating equation;
      V_i<-WCov(A_i,mujk,m_i,L,invT)



      R_hat<-((mY_i-m_i)/sqrt(m_i*(1-m_i)))%*%t((mY_i-m_i)/sqrt(m_i*(1-m_i)))

      H_i1<-D%*%SymmInv(omega)%*%t(D)%*%solve(V_i)
      I_ni<-diag(ncol(H_i1))
      #for(kkk in 1:n[i]){
      #I_element<-matrix(rep(1,L^2),L)
      #  if (kkk==1) {I_ni<-I_element}else {I_ni<-adiag(I_ni,I_element)}
      #}
      #MAEE<-G_list[[i]]%*%(mY_i-m_i)%*%t(mY_i-m_i)

      #####create bias correct matrix (MAEE)
      id_seq<-seq(1:nrow(H_i1))
      id_row<-Star(id_seq,L)
      id_col<-DStar(id_seq,L)

      lll<-NULL
      R_tilde<-NULL
      for (lll in 1:length(id_row)) {

        G_ij<-solve(I_ni-H_i1)[id_row[lll],]
        R_hat_ij<-R_hat[,id_col[lll]]
        R_tilde_ij<-G_ij%*%R_hat_ij
        #R_tilde_ij<-MAEE[id_row[lll],id_col[lll]]
        if (lll==1){R_tilde<-R_tilde_ij}else{R_tilde<-c(R_tilde,R_tilde_ij)}

      }

      ###Orthogonalized residuals;
      b1<-(mujk*(1-muk)*(muk-mujk))/(muj*(1-muj)*muk*(1-muk)-(mujk-muj*muk)*(mujk-muj*muk))
      b2<-(mujk*(1-muj)*(muj-mujk))/(muj*(1-muj)*muk*(1-muk)-(mujk-muj*muk)*(mujk-muj*muk))


      if (MMORTH==TRUE){
        R<-R_tilde*sqrt(muj*(1-muj))*sqrt(muk*(1-muk))-(mujk-muj*muk)-(b1-muk)*(Ystar-muj)-(b2-muj)*(Ydstar-muk)###this R here is MMORTH
      }else{
        R<-(Ystar*Ydstar)-mujk-b1*(Ystar-muj)-b2*(Ydstar-muk)
      }


      ###Va = WACov( mujk, muj, muk, c);
      #VEE<-(muj-mujk)*mujk*(muk-mujk)*(-1+muj+muk-mujk)/(mujk*mujk-2*muj*muk*mujk+muj*muk*(muj+muk-1))

      ###For Cov(R) ;
      ll<-1
      for (j in 1:(n[i]-1)) {
        for (k in (j+1):n[i]) {
          Vr<-matrix(0L,nrow = L^2,ncol = L^2)
          #Vr_j<-matrix(0L,nrow = L^2,ncol = L^2)
          index<-c(((ll-1)*L^2+1):(ll*L^2))
          for (rrow in 1:L^2) {
            for (ccol in 1:rrow) {
              mujk_s<-mujk[index]
              muj_s<-muj[index]
              muk_s<-muk[index]
              #b1_s<-b2[index]###is this correct???
              b1_s<-b1[index]
              b2_s<-b2[index]

              Vr[rrow,ccol]<-Rcov(mujk_s,muj_s,muk_s,b1_s,b2_s,L,rrow,ccol)
              #Vr_j[rrow,ccol]<-Rcov(mujk_s,muj_s,muk_s,b1_s,b2_s,L,rrow,ccol)

              Vr[ccol,rrow]<-Vr[rrow,ccol]
              #Vr_j[ccol,rrow]<-Vr_j[rrow,ccol]
            }
          }
          if (j==1 & k==2){Va<-Vr}else{Va<-magic::adiag(Va,Vr)}

          ll<-ll+1
        }
      }
      ###Da := partial xi / partial alpha = ;
      a<-1/mujk+1/(1-muj-muk+mujk)+1/(muj-mujk)+1/(muk-mujk)
      dmujk_da<-VecDiv(z_i,a)
      Da<-dmujk_da

      ###Va = diag(VEE);
      iVa<-solve(Va)
      U_i[(p+1):(p+q)]<-t(Da)%*%iVa%*%R
      DVD_i[(p+1):(p+q),(p+1):(p+q)]<-t(Da)%*%iVa%*%Da

      ###Working covariance for Beta estimating equation;
      #V_i<-WCov(A_i,mujk,m_i,L,invT)
    }
    ###If data not independent;
    invV_i<-solve(V_i)

    U_b<-t(D)%*%invV_i%*%(y_i-p_i)
    U_i[1:p]<-U_b
    DVD_i[1:p,1:p]<-t(D)%*%invV_i%*%D

    U<-U+U_i
    UUt<-UUt+U_i%*%t(U_i)
    DVD<-DVD+DVD_i

    #print(i)
    #print(DVD_i)
    #print(DVD)
    #record elements for bias correction: D, invV, Da, iVa
    D_list[[i]]<-D
    iV_list[[i]]<-invV_i
    mu_list[[i]]<-(y_i-p_i)
    if (independence==0){
      Da_list[[i]]<-Da
      iVa_list[[i]]<-iVa
      R_list[[i]]<-R
    } else {
      Da_list[[i]]<-NULL
      iVa_list[[i]]<-NULL
      R_list[[i]]<-NULL
    }

  }

  #bias correction
  i<-NULL
  if (independence==0){
    for (i in 1:length(n)){
      U_BC0_i<-rep(0,p+q)
      U_BC1a_i<-rep(0,p+q)
      U_BC1b_i<-rep(0,p+q)
      U_BC2_i<-rep(0,p+q)
      U_BC3_i<-rep(0,p+q)
      #calculate cluster leverage matrix H
      H_1i<-D_list[[i]]%*%solve(DVD[1:p,1:p])%*%t(D_list[[i]])%*%iV_list[[i]]
      H_2i<-Da_list[[i]]%*%solve(DVD[(p+1):(p+q),(p+1):(p+q)])%*%t(Da_list[[i]])%*%iVa_list[[i]]

      I_1i<-diag(ncol(H_1i))
      I_2i<-diag(ncol(H_2i))

      BC0_1i<-I_1i
      BC0_2i<-I_2i

      U_BC0_i[1:p]<-t(D_list[[i]])%*%iV_list[[i]]%*%BC0_1i%*%mu_list[[i]]
      U_BC0_i[(p+1):(p+q)]<-t(Da_list[[i]])%*%iVa_list[[i]]%*%BC0_2i%*%R_list[[i]]

      UU_BC0_i<-U_BC0_i%*%t(U_BC0_i)

      BC2_1i<-solve(I_1i-H_1i)
      BC2_2i<-solve(I_2i-H_2i)

      U_BC1a_i[1:p]<-t(D_list[[i]])%*%iV_list[[i]]%*%BC2_1i%*%mu_list[[i]]
      U_BC1b_i[1:p]<-t(D_list[[i]])%*%iV_list[[i]]%*%mu_list[[i]]

      U_BC1a_i[(p+1):(p+q)]<-t(Da_list[[i]])%*%iVa_list[[i]]%*%BC2_2i%*%R_list[[i]]
      U_BC1b_i[(p+1):(p+q)]<-t(Da_list[[i]])%*%iVa_list[[i]]%*%R_list[[i]]

      UU_BC1_1i<-U_BC1a_i%*%t(U_BC1b_i)
      UU_BC1_2i<-U_BC1b_i%*%t(U_BC1a_i)

      U_BC2_i[1:p]<-t(D_list[[i]])%*%iV_list[[i]]%*%BC2_1i%*%mu_list[[i]]
      U_BC2_i[(p+1):(p+q)]<-t(Da_list[[i]])%*%iVa_list[[i]]%*%BC2_2i%*%R_list[[i]]
      UU_BC2_i<-U_BC2_i%*%t(U_BC2_i)

      C_1ii<-diag(t(D_list[[i]])%*%iV_list[[i]]%*%D_list[[i]]%*%solve(DVD[1:p,1:p]))
      C_2ii<-diag(t(Da_list[[i]])%*%iVa_list[[i]]%*%Da_list[[i]]%*%solve(DVD[(p+1):(p+q),(p+1):(p+q)]))

      for (j1 in 1:length(C_1ii)) {
        C_1ii[j1]<-min(0.75,C_1ii[j1])
      }
      for (j2 in 1:length(C_2ii)) {
        C_2ii[j2]<-min(0.75,C_2ii[j2])
      }

      BC3_1i<-diag((1-C_1ii)^(-0.5))
      BC3_2i<-diag((1-C_2ii)^(-0.5))

      U_BC3_i[1:p]<-BC3_1i%*%t(D_list[[i]])%*%iV_list[[i]]%*%mu_list[[i]]
      U_BC3_i[(p+1):(p+q)]<-BC3_2i%*%t(Da_list[[i]])%*%iVa_list[[i]]%*%R_list[[i]]
      UU_BC3_i<-U_BC3_i%*%t(U_BC3_i)

      if (i==1){
        UU_BC0<-UU_BC0_i
        UU_BC1_1<-UU_BC1_1i
        UU_BC1_2<-UU_BC1_2i
        UU_BC2<-UU_BC2_i
        UU_BC3<-UU_BC3_i
      }else{
        UU_BC0<-UU_BC0+UU_BC0_i
        UU_BC1_1<-UU_BC1_1+UU_BC1_1i
        UU_BC1_2<-UU_BC1_2+UU_BC1_2i
        UU_BC2<-UU_BC2+UU_BC2_i
        UU_BC3<-UU_BC3+UU_BC3_i
      }
    }
  } else {
    for (i in 1:length(n)){
      U_BC0_i<-rep(0,p)
      U_BC1a_i<-rep(0,p)
      U_BC1b_i<-rep(0,p)
      U_BC2_i<-rep(0,p)
      U_BC3_i<-rep(0,p)
      H_1i<-D_list[[i]]%*%solve(DVD[1:p,1:p])%*%t(D_list[[i]])%*%iV_list[[i]]

      I_1i<-diag(ncol(H_1i))

      BC0_1i<-I_1i

      U_BC0_i[1:p]<-t(D_list[[i]])%*%iV_list[[i]]%*%BC0_1i%*%mu_list[[i]]

      UU_BC0_i<-U_BC0_i%*%t(U_BC0_i)

      BC2_1i<-solve(I_1i-H_1i)

      U_BC1a_i[1:p]<-t(D_list[[i]])%*%iV_list[[i]]%*%BC2_1i%*%mu_list[[i]]
      U_BC1b_i[1:p]<-t(D_list[[i]])%*%iV_list[[i]]%*%mu_list[[i]]


      UU_BC1_1i<-U_BC1a_i%*%t(U_BC1b_i)
      UU_BC1_2i<-U_BC1b_i%*%t(U_BC1a_i)

      U_BC2_i[1:p]<-t(D_list[[i]])%*%iV_list[[i]]%*%BC2_1i%*%mu_list[[i]]
      UU_BC2_i<-U_BC2_i%*%t(U_BC2_i)

      C_1ii<-diag(t(D_list[[i]])%*%iV_list[[i]]%*%D_list[[i]]%*%solve(DVD[1:p,1:p]))

      for (j1 in 1:length(C_1ii)) {
        C_1ii[j1]<-min(0.75,C_1ii[j1])
      }


      BC3_1i<-diag((1-C_1ii)^(-0.5))

      U_BC3_i[1:p]<-BC3_1i%*%t(D_list[[i]])%*%iV_list[[i]]%*%mu_list[[i]]
      UU_BC3_i<-U_BC3_i%*%t(U_BC3_i)

      if (i==1){
        UU_BC0<-UU_BC0_i
        UU_BC1_1<-UU_BC1_1i
        UU_BC1_2<-UU_BC1_2i
        UU_BC2<-UU_BC2_i
        UU_BC3<-UU_BC3_i
      }else{
        UU_BC0<-UU_BC0+UU_BC0_i
        UU_BC1_1<-UU_BC1_1+UU_BC1_1i
        UU_BC1_2<-UU_BC1_2+UU_BC1_2i
        UU_BC2<-UU_BC2+UU_BC2_i
        UU_BC3<-UU_BC3+UU_BC3_i
      }
    }
  }

  UU_BC1<-(UU_BC1_1+UU_BC1_2)/2

  Result_score<-list(U,UUt,DVD,UU_BC0,UU_BC1,UU_BC2,UU_BC3)
  return(Result_score)
}
