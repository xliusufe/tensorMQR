
##--------------Estimation without Penalty----------------------##
mqr_cv <- function(Y,X,ncv,r1_index,r3_index,S,U,V,mu,isSym,opts){
  n = opts$n
  p = opts$p
  q = opts$q
  len_cv = ceiling(n/ncv)  
  RSS = rep(0,length(r1_index)*length(r3_index))
  for(jj in 1:ncv){ # start CV
    cv.id = ((jj-1)*len_cv+1):(jj*len_cv)
    if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
    Ytrain = Y[-cv.id,]
    Xtrain = X[-cv.id,]
    Ytest = Y[cv.id,]
    Xtest = X[cv.id,]
    Ztest = produceZ(Xtest)
    nt = nrow(Ytest)
    RSS0 = NULL
    for(r3 in r3_index){
      opts$r3 = r3
      for(r1 in r1_index){
        opts$r1 = r1
        opts$r2 = r1
        if(isSym)
          fit = Estimation(Ytrain,Xtrain,as.matrix(S[1:r3,1:r1^2]),as.matrix(U[,1:r1]),as.matrix(V[,1:r3]),mu,opts)
        else
          fit = EstUnconstr(Y,X,as.matrix(S[1:r3,1:r1^2]),as.matrix(U[,1:r1]),as.matrix(U[,1:r1]),as.matrix(V[,1:r3]),mu,opts)
        if(opts$intercept) RSS0 = c(RSS0,sum((Ytest - matrix(rep(fit$mu,each=nt),nt) - Ztest%*%t(fit$Dnew))^2))
        else RSS0 = c(RSS0,sum((Ytest - Ztest%*%t(fit$Dnew))^2))
      }
    }
    RSS = RSS + RSS0
  } # end of CV
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
