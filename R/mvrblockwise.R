mvrblockwise <- 
  function(Y,X,Z=NULL,method="BIC",ncv=10,penalty="LASSO",isPenColumn=TRUE,group=NULL,lambda=NULL,nlam=50,
           intercept=TRUE,lam_min=1e-6,eps=1e-6,max_step=20,gamma_pen=2,dfmax=NULL,alpha=1){
    
    n <- nrow(Y)
    q <- ncol(Y)
    p <- ncol(X)
    if(is.null(group)) group = rep(1:p)
    gunique <- unique(group)
    G = length(gunique)
    if(G==p){
      eps = eps/sqrt(q)
      fit_mvr <- mvrcolwise(Y,X,Z,method,ncv,penalty,isPenColumn,lambda,nlam,intercept,lam_min,eps,max_step,gamma_pen,dfmax,alpha)
      fit_mvr$group = rep(1,p)
      return(fit_mvr)
    }
    lens = rep(0,G)
    X1 = NULL
    for(g in 1:G){
      lens[g] = sum(group==gunique[g])
      X1 = cbind(X1,X[,which(gunique[g]==group)])
    }
    X2bar = colMeans(X1)
    Ybar = colMeans(Y)
    X2 = X1 - matrix(rep(X2bar,each=n),n)
    Y1 = Y - matrix(rep(Ybar,each=n),n)
    eps = eps/(sqrt(q*max(lens)))
    if(is.null(Z)){
      Zbar = 0
      pz = 0
      Z1 = matrix(0,n,2)
    }
    else{
      Zbar = colMeans(Z)
      pz = ncol(Z)
      Z1 = Z - matrix(rep(Zbar,each=n),n)
      L = solve(chol(t(Z1)%*%Z1/n))
      Z1 = Z1%*%L
    }
    
    if (penalty == "LASSO") pen <- 1
    if (penalty == "MCP")   pen <- 2 
    if (penalty=="SCAD"){    
      gamma_pen <- 3.7
      pen <- 3
    }  
    if (gamma_pen <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
    if (gamma_pen <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
    if (is.null(dfmax)) dfmax = G + 1    
    opts = list(eps=eps,max_step=max_step,n=n,p=G,q=q,pz=pz) 
    opts_pen = list(pen=pen,nlam=nlam,lam_max=1,lam_min=lam_min,gamma_pen=gamma_pen,alpha=alpha,dfmax=dfmax,isPenColumn=isPenColumn) 
    
    if(isPenColumn){
      if (is.null(lambda)) {
        if (nlam < 1) stop("nlambda must be at least 1")
        if (n<=p) lam_min = 1e-2
        setlam = c(1,lam_min,alpha,nlam)
        if(pz){
          fitlm = lm(Y1~Z1-1)
          lambda = setuplambdaMVR_blockwise(fitlm$residuals,X2,nlam,setlam,lens)
        }
        else lambda = setuplambdaMVR_blockwise(Y1,X2,nlam,setlam,lens)
      }
      else  opts_pen$nlam = length(lambda)
    }
    else{
      if (is.null(lambda)) {
        if (nlam < 1) stop("nlambda must be at least 1")
        if (n<=p) lam_min = 1e-2
        setlam = c(1,lam_min,alpha,nlam)
        if(pz){
          fitlm = lm(Y1~Z1-1)
          lambda = setuplambdaMVR_glasso(fitlm$residuals,X2,nlam,setlam,lens)
        }
        else lambda = setuplambdaMVR_glasso(Y1,X2,nlam,setlam,lens)
      }
      else  opts_pen$nlam = nrow(lambda)
    }
    #---------------- The selection by CV or BIC  ---------------------# 
    if(method=="CV") fit_mvr = mvrblockwise_cv(Y1,X2,Z1,ncv,lambda,lens,opts,opts_pen)
    else fit_mvr = mvrblockwise_bic(Y1,X2,Z1,method,lambda,lens,opts,opts_pen)
    if(pz) fit_mvr$Chat = L%*%fit_mvr$Chat
    if(intercept){
      if(pz) fit_mvr$muhat = Ybar-t(fit_mvr$Bhat)%*%X2bar-t(fit_mvr$Chat)%*%Zbar
      else  fit_mvr$muhat = Ybar-t(fit_mvr$Bhat)%*%X2bar
    }
    else fit_mvr$muhat = rep(0,q)
    
    Bhat = fit_mvr$Bhat
    d = 0
    for(g in 1:G){
      Bhat[which(gunique[g]==group),] = fit_mvr$Bhat[(d+1):(d+lens[g]),]
      d = d+lens[g]
    }
    fit_mvr$Bhat = Bhat
    fit_mvr$group = group
    return(fit_mvr)
  }