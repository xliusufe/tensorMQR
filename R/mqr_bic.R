
##--------------main by BIC without sparsity----------------------##
mqr_bic <- function(Y,X,r1_index,r3_index,S,U,V,opts){
  n = opts$n
  p = opts$p
  q = opts$q
  
  # temp = kronecker(U,U)
  # alpha = (Y-Z%*%temp%*%t(V%*%S))%*%(V%*%S)
  # print(2*sum(alpha*(Z%*%temp)))
  
  
  
  RSS = NULL
  for(r3 in r3_index){
    opts$r3 = r3
    for(r1 in r1_index){
      opts$r1 = r1
      opts$r2 = r1
      fit = Estimation(Y,X,as.matrix(S[1:r3,1:r1^2]),as.matrix(U[,1:r1]),as.matrix(V[,1:r3]),opts)
      #df = r1*r1*r3+2*p*r1+q*r3-2*r1^2-r3^2/2
      df = r1*(r1+1)*r3/2+p*r1+q*r3-r1^2-r3^2/2
      RSS = c(RSS,2*log(fit$likhd)+log(n)*df/n)
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
  fit = Estimation(Y,X,as.matrix(S[1:r3_opt,1:r1_opt^2]),as.matrix(U[,1:r1_opt]),as.matrix(V[,1:r3_opt]),opts)  
  return(list(Dnew=fit$Dnew, 
              rss=fit$likhd,
              rk_opt=c(r1_opt,r3_opt),
              selected=selected,
              Y = Y,
              X = X
              )
         )
}