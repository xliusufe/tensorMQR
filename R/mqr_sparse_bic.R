
##--------------Estimation with Penalty by BIC----------------------##
mqr_sparse_bic <- function(Y,X,method,r1_index,r3_index,S,U,V,lambda,isSym,mu,opts,opts_pen){    
  p = opts$p
  q = opts$q
  n = opts$n
  nlam = opts_pen$nlam
  RSS = NULL
  for(r3 in r3_index){
    opts$r3 = r3
    for(r1 in r1_index){
      opts$r1 = r1
      opts$r2 = r1
      if(opts_pen$isPenColumn){
        if(isSym)
          fit = EstPenColumn(Y,X,as.matrix(S[1:r3,1:(r1^2)]),as.matrix(U[,1:r1]),as.matrix(V[,1:r3]),lambda,opts,opts_pen)
        else
          fit = EstUnconstrPen(Y,X,as.matrix(S[1:r3,1:(r1^2)]),as.matrix(U[,1:r1]),as.matrix(U[,1:r1]),
                               as.matrix(V[,1:r3]),lambda,mu,opts,opts_pen)
        df = fit$df*r1
      }
      else{
        fit = EstPenSingle(Y,X,as.matrix(S[1:r3,1:r1^2]),as.matrix(U[,1:r1]),as.matrix(V[,1:r3]),lambda,opts,opts_pen) 
        df = r1*(r1+1)*r3/2 + colSums(fit$betapath) + q*r3 - r1^2-r3^2/2
      }
      loglikelih = (n*q)*log(fit$likhd/(n*q))
      bic <- switch (method,
                     BIC = loglikelih + log(n*q)*df,
                     AIC = loglikelih + 2*df,
                     GCV = fit$likhd*(n*q)/(n*q-df)^2,
                     EBIC = loglikelih + log(n*q)*df + 2*(lgamma(q*p*(p-1)/2+1) 
                                       - lgamma(df+1) - lgamma(q*p*(p-1)/2-df+1))
      )      
      RSS = c(RSS,bic)
      }
    }

  selected = which.min(RSS)
  qj = ceiling(selected/nlam)
  qj1 = selected%%nlam
  if(qj1==0) qj1=nlam
  
  lambda_opt = lambda[qj1]
  opt = assig(c(length(r1_index),length(r3_index)))[,qj]
  r1_opt = r1_index[opt[1]]
  r3_opt = r3_index[opt[2]]
  #---------------- The estimation after selection ---------------------#
  opts$r1 = r1_opt
  opts$r2 = r1_opt
  opts$r3 = r3_opt
  if(opts_pen$isPenColumn){
    if(isSym) fit_opt = EstPenColumn(Y,X,as.matrix(S[1:r3_opt,1:(r1_opt^2)]),as.matrix(U[,1:r1_opt]),as.matrix(V[,1:r3_opt]),lambda,opts,opts_pen)
    else fit_opt = EstUnconstrPen(Y,X,as.matrix(S[1:r3_opt,1:(r1_opt^2)]),as.matrix(U[,1:r1_opt]),as.matrix(U[,1:r1_opt]),as.matrix(V[,1:r3_opt]),lambda,mu,opts,opts_pen)
    activeF = activeX = fit_opt$betapath[,qj1]
  }
  else{
    fit_opt = EstPenSingle(Y,X,as.matrix(S[1:r3_opt,1:(r1_opt^2)]),as.matrix(U[,1:r1_opt]),as.matrix(V[,1:r3_opt]),lambda,opts,opts_pen)
    activeF = matrix(fit_opt$betapath[,qj1],q,p-1)
    activeX = fit_opt$activeXpath[,qj1]
    print(dim(fit_opt$betapath))
    print(length(fit_opt$betapath[,qj1]))
  }
  Unew=matrix(fit_opt$Upath[,qj1],nrow=p)
  Vnew=matrix(fit_opt$Vpath[,qj1],nrow=q)
  Snew=matrix(fit_opt$Spath[,qj1],nrow=r3_opt)  
  return(list(rss=fit_opt$likhd[qj1],
              df = fit_opt$df,
              activeF = activeF,
              activeX = activeX,
              lambda = lambda,
              selectedID = selected,
              lambda_opt=lambda_opt,
              RSS = RSS,
              Unew=Unew,
              Vnew=Vnew,
              Snew=Snew,
              rk_opt=c(r1_opt,r3_opt),
              Y = Y,
              X = X
              )
         )
}
