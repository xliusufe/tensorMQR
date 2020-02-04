# list of functions which are used by R codes.

assig <- function(n_args){
  cargs <- vector("list", length(n_args))
  for(i in 1:length(n_args)) cargs[[i]] <- 1:n_args[i]
  t(expand.grid(cargs))
}

##--------------generate data----------------------##
generateData <- function(n,q,p,s,D3,SigmaX=NULL,SigmaE=NULL,sigma2=NULL,seed_id=NULL){
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s<1) stop("s must be not smaller than 1")
  if(is.null(SigmaX)) SigmaX = diag(p-1)
  if(is.null(SigmaE)) SigmaE = diag(q)
  if(is.null(sigma2)) sigma2 = 0.2
  if(sigma2<=0){
    warning("sigma2 <= 0; set to 0.1")
    sigma2=0.2
  }
  if(is.null(seed_id)) seed_id=1000
  set.seed(seed_id)
  X <- matrix(rnorm(n*(p-1)), nrow = n)%*%chol(SigmaX)
  X <- cbind(rep(1,n),X)
  Z = produceZ(X[,1:s])
  eps <- matrix(rnorm(n*q),n,q)
  Y <- Z%*%t(D3) + sigma2*eps%*%chol(SigmaE)
  return(list(Y=Y,X=X))
}




