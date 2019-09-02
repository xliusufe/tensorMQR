
##--------------Estimation with Penalty by CV----------------------##
mqr_sparse_dr <- 
  function(Y,X,r1_index=NULL,r3_index=NULL,method="BIC",ncv=10,penalty="LASSO",isPenU=0,isPenColumn=1,lambda=NULL,SUV=NULL,
           nlam=50,lam_min=0.001,eps=1e-6,max_step=20,max_step1=20, thresh=1e-4,gamma=2,dfmax=NULL,alpha=1){
    n <- dim(Y)[1]
    q <- dim(Y)[2]
    p <- dim(X)[2]
    if(is.null(r1_index)) r1_index = 1:min(floor(log(n)),p)
    if(is.null(r3_index)) r3_index = 1:min(floor(log(n)),q)
    if (penalty == "LASSO") pen <- 1
    if (penalty == "MCP")   pen <- 2 
    if (penalty=="SCAD"){    
      gamma <- 3
      pen <- 3
    }  
    if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
    if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
    
    if (is.null(dfmax)) dfmax = p + 1
    
    r1_max = max(r1_index)
    r2_max = r1_max
    r3_max = max(r3_index)
    opts = list(utol=1e-4,ftol=eps,Pitol=1e-4,tau_min=1e-3,eta=0.1,tiny=1e-13,gamma=0.85,rhols=1e-4,
                max_step=max_step,max_step1=max_step1,is_LR=1,n=n,r1=r1_max,r2=r2_max,r3=r3_max,p=p,q=q,onStiefel=1) # 18
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
    if (is.null(lambda)) {
      if (nlam < 1) stop("nlambda must be at least 1")
      if (n<=p) lam_min = 1e-2
      setlam = c(1,lam_min,alpha,nlam)
      fit = Estimation(Y,X,S,U,V,opts)
      S = fit$S; U = fit$U; V = fit$V;
      lambda = setuplambda(Y,X,S,U,V,isPenU,nlam,setlam)
    }
    else  nlam = length(lambda)
    
    opts_pen = list(pen=pen,nlam=nlam,lam_max=1,lam_min=lam_min,gamma_pen=gamma,alpha=alpha,dfmax=dfmax,gamma_tanh=1000,
                    thresh=thresh,isPenU=isPenU,isPenColumn=isPenColumn,delta=1,max_step=5,max_step1=10,isFISC=1) #15
    #---------------- The selection by CV  ---------------------#  
    if(method=="BIC") fit_dr = mqr_sparse_bic(Y,X,r1_index,r3_index,S,U,V,lambda,opts,opts_pen)
    if(method=="CV") fit_dr = mqr_sparse_cv(Y,X,ncv,r1_index,r3_index,S,U,V,lambda,opts,opts_pen)
    
    return(fit_dr)
  }