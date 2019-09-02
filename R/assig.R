# list of functions which are used by R codes.

assig <- function(n_args){
  cargs <- vector("list", length(n_args))
  for(i in 1:length(n_args)) cargs[[i]] <- 1:n_args[i]
  t(expand.grid(cargs))
}

##--------------generate data----------------------##
generateData <- function(n,q,s,p,D3,SigmaX=diag(p-1),sigma2=0.2,seed_id=1000, rho=0.0){
  set.seed(seed_id)
  X <- matrix(rnorm(n*(p-1)), nrow = n)%*%chol(SigmaX)
  X <- cbind(rep(1,n),X)
  Z = produceZ(X[,1:s])
  Sigma = matrix(rho,q,q)+diag(1-rho,q,q)
  A = chol(Sigma) 
  eps <- matrix(rnorm(n*q),n,q)
  Y <- Z%*%t(D3) + sigma2*eps%*%A
  return(list(Y=Y,X=X))
}




