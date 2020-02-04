
##--------------Estimation with Penalty by CV----------------------##
mqr_sparse_dr <- 
  function(Y,X0,r1_index=NULL,r3_index=NULL,method="BIC",ncv=10,penalty="LASSO",isPenU=0,isPenColumn=TRUE,
           lambda=NULL,SUV=NULL,isSym=TRUE,initMethod="LASSO",nlam=50,lam_min=0.001,ftol=1e-6,max_step=20,
           max_step1=20,eps=1e-4,thresh=1e-4,gamma_pen=2,dfmax=NULL,alpha=1){
    n <- nrow(Y)
    q <- ncol(Y)

    if(is.null(initMethod)){
        p = ncol(X0)
        Z = produceZ(X0)
        selectX = rep(1,p)
        fit_mvr = NULL
    }
    else{
      p0 <- ncol(X0)
      X2 <- produceX2(X0)[,-1]
      fit_mvr <- mvrcolwise(Y,X2,method="GCV",penalty=initMethod, isPenColumn=TRUE)
      activeX = c(1,fit_mvr$activeX)
      selectX <- selectedP2X(activeX, p0)
      X = X0[,which(selectX==1)]
      p <- ncol(X)
      Z <- produceZ(X)
      SUV <- NULL
    }
    
    if(is.null(r1_index)) r1_index = 1:min(floor(log(n)),p)
    if(is.null(r3_index)) r3_index = 1:min(floor(log(n)),q)
    if (penalty == "LASSO") pen <- 1
    if (penalty == "MCP")   pen <- 2 
    if (penalty=="SCAD"){    
      gamma_pen <- 3.7
      pen <- 3
    }  
    if (gamma_pen <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
    if (gamma_pen <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
    
    if (is.null(dfmax)) dfmax = p + 1
    r1_max = max(r1_index)
    r2_max = r1_max
    r3_max = max(r3_index)
    opts = list(eps=eps,eps1=eps,utol=1e-4,ftol=ftol,Pitol=1e-4,tau_min=1e-3,eta=0.1,tiny=1e-13,gamma=0.85,rhols=1e-4,
                max_step=max_step,max_step1=max_step1,isLR=1,n=n,r1=r1_max,r2=r2_max,r3=r3_max,p=p,q=q,onStiefel=1) 
    # initial U,V,S
    if(is.null(SUV)){
      set.seed(1)
      U = rbind(diag(r1_max), matrix(0,p-r1_max,r1_max))
      V = rbind(diag(r3_max), matrix(0,q-r3_max,r3_max))
      S = matrix(rnorm(r1_max*r2_max*r3_max),r3_max,r1_max*r2_max)
    }
    else{
      U = SUV$U 
      V = SUV$V
      S = SUV$S
    } 
    fit = Estimation(Y,Z,S,U,V,opts)
    S1 = fit$S; U1 = fit$U; V1 = fit$V
    if (is.null(lambda)) {
      if (nlam < 1) stop("nlambda must be at least 1")
      if (n<=p) lam_min = 1e-2
      setlam = c(1,lam_min,alpha,nlam)
      lambda = setuplambda(Y,Z,S1,U1,V1,isPenU,nlam,setlam)
    }
    else  nlam = length(lambda)
    opts_pen = list(pen=pen,nlam=nlam,lam_max=1,lam_min=lam_min,gamma_pen=gamma_pen,alpha=alpha,dfmax=dfmax,gamma_tanh=1000,
                    thresh=thresh,isPenU=isPenU,isPenColumn=isPenColumn,delta=1,max_step=5,max_step1=10,isFISC=1) #15
    #---------------- The selection by CV  ---------------------#  
    if(method=="CV") fit_dr = mqr_sparse_cv(Y,Z,ncv,r1_index,r3_index,S1,U1,V1,lambda,isSym,opts,opts_pen)
    else fit_dr = mqr_sparse_bic(Y,Z,method,r1_index,r3_index,S1,U1,V1,lambda,isSym,opts,opts_pen)
    
    fit_dr$Y = Y
    fit_dr$X = X0
    fit_dr$Xselect=which(selectX==1)
    return(fit_dr)
  }