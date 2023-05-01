
#' function: ORTH.Ord
#'
#' This function is designed for analyzing correlated ordinal data with ability to correct small-sample bias.
#'
#' The method a modified version of alternating logistic regressions with estimation based on orthogonalized residuals (ORTH).
#'   The within-cluster association between ordinal responses is modeled by global pairwise odds ratios (POR).
#'   A small-sample bias correction to POR parameter estimates based on matrix multiplicative adjusted orthogonalized residuals (MMORTH)
#'   for correcting estimating equations, and bias-corrected sandwich estimators with different options for covariance estimation,
#'   i.e. BC1 (Kauermann & Zeger (1986)), BC2 (Mancl & DeRouen (2001)), and BC3 (Fay & Graubard (2001)).
#'
#' @param formula_mean the symbolic description of the marginal mean model that contains the ordinal outcome and marginal mean covariates.
#' @param data_mean the data set containing the ordinal outcome and marginal mean covariates.
#' @param cluster cluster ID (consecutive integers) in data_mean.
#' @param formula_por the symbolic description of marginal association model in the form of a one-sided formula, default is NULL. When leaving formula_por as default, independence working correlation will be used.
#' @param data_por a data set for marginal association model, default is NULL. When leaving data_por as default, independence working correlation will be used.
#' @param MMORTH a logical value to indicate if matrix-adjusted estimating equations will be applied for the association estimation, default is FALSE.
#' @param BC an option to apply bias-correction on covariance estimation, default is NULL. Possible values are "BC1", "BC2", or "BC3".
#' @param init_beta pre-specified starting values for parameters in the mean model, default is NULL.
#' @param init_alpha pre-specified starting values for parameters in the association model, default is NULL.
#' @param miter maximum number of iterations for Fisher scoring, default is 30.
#' @param crit_level tolerance for convergence, default is 0.0001.
#' @return a matrix for estimation results
#' @export


ORTH.Ord<-function(formula_mean,data_mean,cluster,formula_por=NULL,data_por=NULL,MMORTH=FALSE,BC=NULL,init_beta=NULL,init_alpha=NULL,miter=30,crit_level=0.0001){
  #require(MASS)
  #require(matrixcalc)
  #require(magic)

  if (missing(formula_mean)) { stop("You need to provide a model formula for the mean.")        }
  if (missing(data_mean))    { stop("You need to provide data for marginal mean effects.")      }
  if (missing(cluster))      { stop("You need to provide the cluster id.") }
  if (length(formula_mean)!=3){stop("Formula for marginal mean model must be two-sided ")} #length of two-sided formula is 3

  formula1<-stats::as.formula(formula_mean)
  name_mean<-all.vars(formula1)
  response<-name_mean[1]
  covar_mean<-name_mean[2:length(name_mean)]

  #check if all variables for the models are numeric
  numeric_chk<-sapply(data_mean[,c(paste0(response),paste0(covar_mean))], is.numeric)
  if(all(numeric_chk)==F){stop("The variables specified in marginal mean model must be all numeric.")}

  #if initial values of beta are specified, the length must be consistent with the number of covariates in mean model
  if (length(init_beta)>0){
    if (length(init_beta)!=(length(name_mean)-1)){stop("The number of initial values for beta must be consistent with the number of covariates in the mean model.")}
  }

  #check number of observations and remove missing data points
  row_n1<-dim(data_mean)[1]
  data_mean<-data_mean[stats::complete.cases(data_mean[,c(response,covar_mean)]),]
  row_n2<-dim(data_mean)[1]
  cat("The data for the marginal mean model \n")
  cat("Number of observations read:",row_n1,"\n")
  cat("Number of observations used:",row_n2,"\n")

  ordinal_L<-length(levels(as.factor(c(data_mean[,paste0(response)]))))-1#check how many levels the ordinal response has, must be >=2
  if (ordinal_L>=2){
    cat("The ordinal response has ",ordinal_L+1," levels. \n")}else{
      stop("The ordinal response must have 3 levels or more")
    }

  #check number of cluster and number of rows for correlation data
  id_vector<-data_mean[,cluster]
  id_dim<-table(id_vector)
  cluster_n<-dim(id_dim)
  cat("Number of clusters is:",cluster_n,"\n")
  for (i in 1:dim(id_dim)) {
    row_n3<-choose(id_dim[[i]],2)
    if (i==1){row_n4<-row_n3}else{row_n4<-row_n4+row_n3}
  }
  #row_n4 is the number of pairwise observations expected in data_por

  if (length(data_por)==0 | length(formula_por)==0){
    independence<-1
    covar_por<-NULL
    init_alpha<-NULL#if independence working correlation is used, no pre-specified initial values
  } else {
    independence<-0
    if (length(formula_por)!=2){stop("Formula for marginal correlation model must be one-sided.")}
    if (dim(data_por)[1]!=row_n4){
      cat("The number of observations in data_por is:",dim(data_por)[1],"\n")
      cat("However, based on the number of clusters and observations per cluster, \n")
      cat("the number of observations in data_por should be:",row_n4,"\n")
      stop("The number of observations in data_por is incorrect, please check your data for both mean model and correlation model.")}
    formula2<-stats::as.formula(formula_por)
    covar_por<-all.vars(formula2)
    #check if all variables for the models are numeric
    numeric_chk<-sapply(data_por[,c(paste0(covar_por))], is.numeric)
    if(all(numeric_chk)==F){stop("The variables specified in the correlation model must be all numeric.")}
    #if initial values of alpha are specified, the length must be consistent with the number of covariates in correlation model
    if (length(init_alpha)>0){
      if (length(init_alpha)!=length(covar_por)){stop("The number of initial values for alpha must be consistent with the number of covariates in the correlation model.")}
    }
  }

  data<-MakeData(response,covar_mean,covar_por,data_mean,data_por,cluster,ordinal_L,independence)

  results<-Fit(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],data[[6]],covar_por,independence,MMORTH,init_beta,init_alpha,miter,crit_level)

  if(all(diag(results[[11]])>=0) & all(diag(results[[12]])>=0) & all(diag(results[[13]])>=0) &all(diag(results[[14]])>=0)){
    p<-length(results[[1]])
    q<-length(results[[2]])
    result_cov<-results[[3]]
    result_beta<-results[[1]]
    result_alpha<-results[[2]]

    #BC
    result_cov_BC0<-results[[11]]
    result_cov_BC1<-results[[12]]
    result_cov_BC2<-results[[13]]
    result_cov_BC3<-results[[14]]


    result_stderr<-sqrt(diag(result_cov[1:p,1:p]))

    #result_chisq<-(result_beta/result_stderr)^2
    #result_pvalue<-1-pchisq(result_chisq,1)

    result_Z<-(result_beta/result_stderr)
    result_pvalue<-(1-stats::pnorm(abs(result_Z)))*2

    if (independence==0){
      result_stderr_a<-sqrt(diag(result_cov[(p+1):(p+q),(p+1):(p+q)]))

      #result_chisq_a<-(result_alpha/result_stderr_a)^2
      #result_pvalue_a<-1-pchisq(result_chisq_a,1)

      result_Z_a<-(result_alpha/result_stderr_a)
      result_pvalue_a<-(1-stats::pnorm(abs(result_Z_a)))*2
    }


    #BC0
    #result_stderr_BC0<-sqrt(diag(result_cov_BC0[1:p,1:p]))


    #if (independence==0){
    #  result_stderr_a_BC0<-sqrt(diag(result_cov_BC0[(p+1):(p+q),(p+1):(p+q)]))
    #}

    #BC1
    result_stderr_BC1<-sqrt(diag(result_cov_BC1[1:p,1:p]))

    result_Z_BC1<-(result_beta/result_stderr_BC1)
    result_pvalue_BC1<-(1-stats::pnorm(abs(result_Z_BC1)))*2


    if (independence==0){
      result_stderr_a_BC1<-sqrt(diag(result_cov_BC1[(p+1):(p+q),(p+1):(p+q)]))

      result_Z_a_BC1<-(result_alpha/result_stderr_a_BC1)
      result_pvalue_a_BC1<-(1-stats::pnorm(abs(result_Z_a_BC1)))*2

    }

    #BC2
    result_stderr_BC2<-sqrt(diag(result_cov_BC2[1:p,1:p]))

    result_Z_BC2<-(result_beta/result_stderr_BC2)
    result_pvalue_BC2<-(1-stats::pnorm(abs(result_Z_BC2)))*2


    if (independence==0){
      result_stderr_a_BC2<-sqrt(diag(result_cov_BC2[(p+1):(p+q),(p+1):(p+q)]))

      result_Z_a_BC2<-(result_alpha/result_stderr_a_BC2)
      result_pvalue_a_BC2<-(1-stats::pnorm(abs(result_Z_a_BC2)))*2

    }

    #BC3
    result_stderr_BC3<-sqrt(diag(result_cov_BC3[1:p,1:p]))

    result_Z_BC3<-(result_beta/result_stderr_BC3)
    result_pvalue_BC3<-(1-stats::pnorm(abs(result_Z_BC3)))*2


    if (independence==0){
      result_stderr_a_BC3<-sqrt(diag(result_cov_BC3[(p+1):(p+q),(p+1):(p+q)]))

      result_Z_a_BC3<-(result_alpha/result_stderr_a_BC3)
      result_pvalue_a_BC3<-(1-stats::pnorm(abs(result_Z_a_BC3)))*2

    }



    final_beta_est<-data.frame(rbind(result_beta,result_stderr,result_Z,result_pvalue,result_stderr_BC1,result_Z_BC1,result_pvalue_BC1,result_stderr_BC2,result_Z_BC2,result_pvalue_BC2,result_stderr_BC3,result_Z_BC3,result_pvalue_BC3))

    for (i in 1:results[[6]]){
      intcpt<-paste0("intercept ",i)
      if (i==1){
        intercept<-intcpt
      } else{
        intercept<-c(intercept,intcpt)
      }
    }

    colnames(final_beta_est)<-c(paste0(intercept),paste0(covar_mean))


    final_beta_est2<-t(final_beta_est)
    colnames(final_beta_est2)<-NULL

    if (independence==0){
      final_alpha_est<-data.frame(rbind(result_alpha,result_stderr_a,result_Z_a,result_pvalue_a,result_stderr_a_BC1,result_Z_a_BC1,result_pvalue_a_BC1,result_stderr_a_BC2,result_Z_a_BC2,result_pvalue_a_BC2,result_stderr_a_BC3,result_Z_a_BC3,result_pvalue_a_BC3))
      colnames(final_alpha_est)<-c(paste0(covar_por))

      final_alpha_est2<-t(final_alpha_est)
      colnames(final_alpha_est2)<-NULL

      final_est<-data.frame(rbind(final_beta_est2,final_alpha_est2))
      names_param<-c(paste0(intercept),paste0(covar_mean),paste0(covar_por))
      final_est<-cbind(names_param,final_est)
    } else {
      final_est<-data.frame(final_beta_est2)
      names_param<-c(paste0(intercept),paste0(covar_mean))
      final_est<-cbind(names_param,final_est)
    }

    rownames(final_est)<-NULL
    colnames(final_est)<-c("Coefficients","Estimate","Std.Error","Z.Value","P.Value","Std.Error(BC1)","Z.Value(BC1)","P.Value(BC1)","Std.Error(BC2)","Z.Value(BC2)","P.Value(BC2)","Std.Error(BC3)","Z.Value(BC3)","P.Value(BC3)")


  }

  if ("BC1" %in% toupper(BC)) {
    BC1<-c("Std.Error(BC1)","Z.Value(BC1)","P.Value(BC1)")
  }else{BC1<-NULL}
  if ("BC2" %in% toupper(BC)) {
    BC2<-c("Std.Error(BC2)","Z.Value(BC2)","P.Value(BC2)")
  }else{BC2<-NULL}
  if ("BC3" %in% toupper(BC)) {
    BC3<-c("Std.Error(BC3)","Z.Value(BC3)","P.Value(BC3)")
  }else{BC3<-NULL}


    final_est2<-final_est[,c("Coefficients","Estimate","Std.Error","Z.Value","P.Value",BC1,BC2,BC3)]

    if (MMORTH==T){cat("MMORTH is used to adjust small-sample bias. \n")}

  return(final_est2)
}


