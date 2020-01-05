mqr_sparse <- 
  function(Y,X,r1=NULL,r3=NULL,method="BIC",ncv=10,penalty="LASSO",isPenU=0,isPenColumn=1,lambda=NULL,SUV=NULL,
           nlam=50,lam_min=1e-3,ftol=1e-6,max_step=20,max_step1=20,eps=1e-4,thresh=1e-4,gamma_pen=2,dfmax=NULL,alpha=1){
    n <- dim(Y)[1]
    q <- dim(Y)[2]
    p <- dim(X)[2]
    isSym=TRUE
    intercept=FALSE
    mu=NULL
    if(is.null(r1)) r1 <- 2 
    if(is.null(r3)) r3 <- 2
    if(r3>q){ 
      r3=q
      warning("r3 should be not greater than q! Reset r3=q.")
    }
    r2 <- r1
    if (penalty == "LASSO") pen <- 1
    if (penalty == "MCP")   pen <- 2 
    if (penalty=="SCAD"){    
      gamma_pen <- 3
      pen <- 3;
    }  
    if (gamma_pen <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
    if (gamma_pen <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
    if (is.null(dfmax)) dfmax = p + 1
    if(!intercept | is.null(mu)) mu = rep(0,q)
    opts = list(eps=eps,eps1=eps,utol=1e-4,ftol=ftol,Pitol=1e-4,tau_min=1e-3,eta=0.1,tiny=1e-13,gamma=0.85,rhols=1e-4,
                max_step=max_step,max_step1=max_step1,isLR=1,n=n,r1=r1,r2=r2,r3=r3,p=p,q=q,onStiefel=1,intercept=intercept) 
    # initial A,B,C,S
    if(is.null(SUV)){
      set.seed(1)
      U = rbind(diag(r1), matrix(0,p-r1,r1))
      V = rbind(diag(r3), matrix(0,q-r3,r3))
      S = matrix(rnorm(r1*r2*r3),r3,r1*r2)
    }
    else{
      U = SUV$U 
      V = SUV$V
      S = SUV$S
    }
    if (is.null(lambda)) {
      if (nlam < 1) stop("nlambda must be at least 1")
      if (n<=p) lam_min = 1e-2
      setlam = c(1,lam_min,alpha,nlam)
      fit = Estimation(Y,X,S,U,V,mu,opts)
      S = fit$S; U = fit$U; V = fit$V;
      lambda = setuplambda(Y,X,S,U,V,isPenU,nlam,setlam)
    }
    else  nlam = length(lambda)
    opts_pen = list(pen=pen,nlam=nlam,lam_max=1,lam_min=lam_min,gamma_pen=gamma_pen,alpha=alpha,dfmax=dfmax,gamma_tanh=1000,
                    thresh=thresh,isPenU=isPenU,isPenColumn=isPenColumn,delta=1,max_step=5,max_step1=10,isFISC=1) #15
    #---------------- The selection by BIC or CV  ---------------------# 
    if(method!="CV"){
      if(isPenColumn)
        if(isSym) fit = EstPenColumn(Y,X,S,U,V,lambda,opts,opts_pen) 
        else fit = EstUnconstrPen(Y,X,S,U,U,V,lambda,mu,opts,opts_pen) 
      else
        fit = EstPenSingle(Y,X,S,U,V,lambda,opts,opts_pen) 
      df = fit$df*r1
      loglikelih =  -n*q * (log(2*pi) + log(fit$likhd))
      bic <- switch (method,
                     BIC = loglikelih + log(n*q)*df,
                     AIC = loglikelih + 2*df,
                     GCV = loglikelih/(1-df/n)^2,
                     EBIC = loglikelih + log(n*q)*df + 2*(lgamma(q*p*(p-1)/2+1) 
                                       - lgamma(df+1) - lgamma(q*p*(p-1)/2-df+1))
      )  
      selected = which.min(bic)
      lambda_opt = lambda[selected]
      activeF = activeX = fit$betapath[,selected]
      Unew=matrix(fit$Upath[,selected],nrow=p)
      Vnew=matrix(fit$Vpath[,selected],nrow=q)
      Snew=matrix(fit$Spath[,selected],nrow=r3)
    }
    if(method=="CV"&&nlam>1){
      len_cv = ceiling(n/ncv)
      RSS = rep(0,nlam)
      for(jj in 1:ncv){
        cv.id = ((jj-1)*len_cv+1):(jj*len_cv)
        if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
        Ytrain = Y[-cv.id,]
        Xtrain = X[-cv.id,]
        Ytest = Y[cv.id,]
        Xtest = X[cv.id,]
        if(isPenColumn){
          if(isSym) fit = EstPenColumn(Ytrain,Xtrain,S,U,V,lambda,opts,opts_pen)
          else  fit = EstUnconstrPen(Ytrain,Xtrain,S,U,U,V,lambda,mu,opts,opts_pen) 
        }
        else
          fit = EstPenSingle(Ytrain,Xtrain,S,U,V,lambda,opts,opts_pen) 
        Dnew = fit$V %*% fit$S %*% t(kroneckerProduct(fit$U, fit$U))
        RSS = RSS + sum((Ytest-produceZ(Xtest)%*%t(Dnew))^2)/2
      } 
      selected = which.min(RSS)
      lambda_opt = lambda[selected]
      
      if(isPenColumn){
        if(isSym)  fit_opt = EstPenColumn(Y,X,as.matrix(S),as.matrix(U),as.matrix(V),lambda[1:selected],opts,opts_pen)
        else  fit_opt = EstUnconstrPen(Y,X,as.matrix(S),as.matrix(U),as.matrix(U),as.matrix(V),lambda[1:selected],mu,opts,opts_pen)
        activeF = activeX = fit_opt$betapath[,selected]
      }
      else{
        fit_opt = EstPenSingle(Y,X,as.matrix(S),as.matrix(U),as.matrix(V),lambda[1:selected],opts,opts_pen)
        activeF = matrix(fit_opt$betapath[,selected],q,p)
        activeX = fit_opt$activeXpath[,selected]
      }
      Unew=matrix(fit_opt$Upath[,selected],nrow=p)
      Vnew=matrix(fit_opt$Vpath[,selected],nrow=q)
      Snew=matrix(fit_opt$Spath[,selected],nrow=r3)
    }
    if(method=="CV"&&nlam==1){
      if(isPenColumn){
        if(isSym)  fit = EstPenColumn(Y,X,S,U,V,lambda,opts,opts_pen)
        else   fit = EstUnconstrPen(Y,X,S,U,U,V,lambda,mu,opts,opts_pen) 
      }
      else
        fit = EstPenSingle(Y,X,S,U,V,lambda,opts,opts_pen) 
      selected = 1
      Unew=matrix(fit$Upath,nrow=p)
      Vnew=matrix(fit$Vpath,nrow=q)
      Snew=matrix(fit$Spath,nrow=r3)
      lambda_opt = lambda
      activeX = activeF = fit$betapath
    }
    return(list(betapath=fit$betapath, 
                rss=fit$likhd[selected],
                df = fit$df,
                lambda = fit$lambda,
                lambda_opt=lambda_opt,
                selectedID = selected,
                activeF = activeF,
                activeX = activeX,
                Unew=Unew,
                Vnew=Vnew,
                Snew=Snew,
                Y = Y,
                X = X
                )
           )
  }