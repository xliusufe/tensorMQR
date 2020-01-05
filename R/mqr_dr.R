
##--------------main by BIC without sparsity----------------------##
mqr_dr <- function(Y,X,method="BIC",ncv=10,r1_index=NULL,r3_index=NULL,SUV=NULL,
                   isSym=TRUE,eps=1e-6,max_step=20,max_step1=20){

  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  intercept=FALSE
  mu=NULL
  if(is.null(r1_index)) r1_index = 1:min(ceiling(log(n)),p)
  if(is.null(r3_index)) r3_index = 1:min(ceiling(log(n)),q)
  #---------------- The selection by BIC  ---------------------#  
  if(is.null(SUV)){
    set.seed(1)
    r1_max = max(r1_index) 
    r3_max = max(r3_index) 
    U = rbind(diag(r1_max), matrix(0,p-r1_max,r1_max))
    V = rbind(diag(r3_max), matrix(0,q-r3_max,r3_max))
    S = matrix(rnorm(r1_max^2*r3_max),r3_max,r1_max^2)
  }
  else{
    U = SUV$U
    V = SUV$V
    S = SUV$S
  }
  if(!intercept | is.null(mu)) mu = rep(0,q)
  if((max(r1_index)>dim(U)[2])|(max(r3_index)>dim(V)[2])){
    r1_index = min(r1_index):min(dim(U)[2],max(r1_index))
    r3_index = min(r3_index):min(dim(V)[2],max(r3_index))
    warning("maximum number of index sequences of r1 and r3 must not be larger than columns of U and V, respectively !")
  }
  opts = list(eps=eps,eps1=eps,utol=1e-4,ftol=eps,Pitol=1e-4,tau_min=1e-3,eta=0.1,tiny=1e-13,gamma=0.85,rhols=1e-4,
              max_step=max_step,max_step1=max_step1,isLR=1,n=n,r1=2,r2=2,r3=2,p=p,q=q,intercept=intercept)
  if(method=="CV")  fit_dr = mqr_cv(Y,X,ncv,r1_index,r3_index,S,U,V,mu,isSym,opts)
  else fit_dr = mqr_bic(Y,X,method,r1_index,r3_index,S,U,V,mu,isSym,opts)
  
  return(fit_dr)
}