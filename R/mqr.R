
##--------------main by BIC without sparsity----------------------##
mqr <- function(Y,X,r1=NULL,r3=NULL,SUV,eps=1e-6,max_step=20){
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  if(is.null(r1)) r1 = 2
  if(is.null(r3)) r3 = 2
  if(is.null(SUV)){
    set.seed(1)
    U = rbind(diag(r1), matrix(0,p-r1,r1))
    V = rbind(diag(r3), matrix(0,q-r3,r3))
    S = matrix(rnorm(r1^2*r3),r3,r1^2)
  }
  else{
    U = SUV$U
    V = SUV$V
    S = SUV$S
  }
  opts = list(utol=1e-4,ftol=eps,Pitol=1e-4,tau_min=1e-3,eta=0.1,tiny=1e-13,gamma=0.85,rhols=1e-4,
              max_step=max_step,max_step1=max_step,is_LR=1,n=n,r1=r1,r2=r1,r3=r3,p=p,q=q)  
  fit = Estimation(Y,X,as.matrix(S),as.matrix(U),as.matrix(V),opts)
  return(list(Dnew=fit$Dnew, 
              rss=fit$likhd,
              S=fit$S,
              U=fit$U,
              V=fit$V,
              Y = Y,
              X = X
  )
  )
}