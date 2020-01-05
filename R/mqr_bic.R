
##--------------main by BIC without sparsity----------------------##
mqr_bic <- function(Y,X,method,r1_index,r3_index,S,U,V,mu,isSym,opts){
  n = opts$n
  p = opts$p
  q = opts$q
  RSS = NULL
  for(r3 in r3_index){
    opts$r3 = r3
    for(r1 in r1_index){
      opts$r1 = r1
      opts$r2 = r1
      if(isSym)
        fit = Estimation(Y,X,as.matrix(S[1:r3,1:r1^2]),as.matrix(U[,1:r1]),as.matrix(V[,1:r3]),mu,opts)
      else
        fit = EstUnconstr(Y,X,as.matrix(S[1:r3,1:r1^2]),as.matrix(U[,1:r1]),as.matrix(U[,1:r1]),as.matrix(V[,1:r3]),mu,opts)
      df = r1*(r1+1)*r3/2+p*r1+q*r3-r1^2-r3^2/2
      loglikelih =  -n*q * (log(2*pi) + log(fit$likhd))
      bic <- switch (method,
                     BIC = loglikelih + log(n*q)*df,
                     AIC = loglikelih + 2*df,
                     GCV = loglikelih/(1-df/n)^2,
                     EBIC = loglikelih + log(n*q)*df + 2*(lgamma(q*p*(p+1)/2+1) 
                                       - lgamma(df+1) - lgamma(q*p*(p+1)/2-df+1))
      )
      RSS = c(RSS,bic)
    }
  }
  selected = which.min(RSS)
  opt = assig(c(length(r1_index),length(r3_index)))[,selected]
  r1_opt = r1_index[opt[1]]
  r3_opt = r3_index[opt[2]]
  #---------------- The estimation after selection ---------------------#
  opts$r1 = r1_opt
  opts$r2 = r1_opt
  opts$r3 = r3_opt
  if(isSym)
    fit = Estimation(Y,X,as.matrix(S[1:r3_opt,1:r1_opt^2]),as.matrix(U[,1:r1_opt]),as.matrix(V[,1:r3_opt]),mu,opts)
  else
    fit = EstUnconstr(Y,X,as.matrix(S[1:r3_opt,1:r1_opt^2]),as.matrix(U[,1:r1_opt]),as.matrix(U[,1:r1_opt]),as.matrix(V[,1:r3_opt]),mu,opts)
  return(list(Dnew=fit$Dnew, 
              rss=fit$likhd,
              mu = fit$mu,
              rk_opt=c(r1_opt,r3_opt),
              selected=selected,
              Y = Y,
              X = X
              )
         )
}