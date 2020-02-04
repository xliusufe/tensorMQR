
##--------------Estimation with Penalty by CV----------------------##
mqr_sparse_cv <- function(Y,X,ncv,r1_index,r3_index,S,U,V,lambda,isSym,opts,opts_pen){
    p <- opts$p
    q <- opts$q
    n <- opts$n
    nlam <- opts_pen$nlam
    len_cv = floor(n/ncv)
    RSS = matrix(0,nlam,length(r1_index)*length(r3_index))
    likhd = rep(0,nlam)
    for(jj in 1:ncv){ # start CV
      cv.id = ((jj-1)*len_cv+1):(jj*len_cv)
      if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
      Ytrain = Y[-cv.id,]
      Xtrain = X[-cv.id,]
      Ytest = Y[cv.id,]
      Xtest = X[cv.id,]
      RSS0 = NULL
      for(r3 in r3_index){
        opts$r3 = r3
        for(r1 in r1_index){
          opts$r1 = r1
          opts$r2 = r1
          if(opts_pen$isPenColumn){
            if(isSym)
              fit = EstPenColumn(Ytrain,Xtrain,as.matrix(S[1:r3,1:(r1^2)]),as.matrix(U[,1:r1]),as.matrix(V[,1:r3]),lambda,opts,opts_pen) 
            else
              fit = EstUnconstrPen(Ytrain,Xtrain,as.matrix(S[1:r3,1:(r1^2)]),as.matrix(U[,1:r1]),as.matrix(U[,1:r1]),
                                   as.matrix(V[,1:r3]),lambda,opts,opts_pen)
            df = r1*(r1+1)*r3/2 + fit$df*r1 + q*r3 - r1^2-r3^2/2
            }
          else{
            fit = EstPenSingleCV(Ytrain,Xtrain,as.matrix(S[1:r3,1:(r1^2)]),as.matrix(U[,1:r1]),as.matrix(V[,1:r3]),lambda,opts,opts_pen)
            df = r1*(r1+1)*r3/2 + colSums(fit$betapath) + q*r3 - r1^2-r3^2/2
          }
          for(kk in 1:nlam){
            Unew=matrix(fit$Upath[,kk],nrow=p)
            Vnew=matrix(fit$Vpath[,kk],nrow=q)
            Snew=matrix(fit$Spath[,kk],nrow=r3) 
            Dnew = Vnew %*% Snew %*% t(kronecker(Unew, Unew))
            likhd[kk] = sum((Ytest-Xtest%*%t(Dnew))^2)/2
          }
          RSS0 = cbind(RSS0,likhd)
        }
      }

      RSS = RSS + RSS0
    } # end of CV
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
      if(isSym)
        fit_opt = EstPenColumn(Y,X,as.matrix(S[1:r3_opt,1:(r1_opt^2)]),as.matrix(U[,1:r1_opt]),as.matrix(V[,1:r3_opt]),lambda,opts,opts_pen) 
      else
        fit_opt = EstUnconstrPen(Y,X,as.matrix(S[1:r3_opt,1:(r1_opt^2)]),as.matrix(U[,1:r1_opt]),as.matrix(U[,1:r1_opt]),
                                 as.matrix(V[,1:r3_opt]),lambda,opts,opts_pen)
      activeF = activeX = fit_opt$betapath[,qj1]
    }
    else{
      fit_opt = EstPenSingle(Y,X,as.matrix(S[1:r3_opt,1:(r1_opt^2)]),as.matrix(U[,1:r1_opt]),as.matrix(V[,1:r3_opt]),lambda,opts,opts_pen)
      activeF = matrix(fit_opt$betapath[,qj1],q,p)
      activeX = fit_opt$activeXpath[,qj1]
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
                RSS = RSS,
                lambda_opt=lambda_opt,
                Unew=Unew,
                Vnew=Vnew,
                Snew=Snew,
                rk_opt=c(r1_opt,r3_opt)
    )
    )
  }