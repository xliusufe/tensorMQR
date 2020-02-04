# tensorMQR1
 Symmetric Tensor Estimation for Quadratic Regression.
 
  For a high-dimensional Multiresponse Quadratic Regression (MQR) with or without aparsity assumptions, 
  treating the coefficients as a third-order symmetric tensor and borrowing Tucker decomposition to reduce the number of parameters.  
  The multivariate sparse group lasso (mcp or scad) and the steepest gradient descent algorithm are used to estimate the tensor for sparsity situation.
# Installation

    #install Rtools 3.5 (http://cran.r-project.org/bin/windows/Rtools)
    #install.packages("devtools")
    #install.packages("Rcpp")
    library(devtools)
    install_github("xliusufe/tensorMQR1")

# Usage

   - [x] [tensorMQR-manual](https://github.com/xliusufe/tensorMQR1/blob/master/inst/tensorMQR1-manual.pdf) ------------ Details of the usage of the package.
# Example

    library(tensorMQR1)

    # Example 1
    # The usage of function "mqr()"
	n <- 200
    p <- 6
    q <- 3	
    D3 <- matrix(runif(q*p^2, 0.7, 1), q, p^2)
    mydata <- generateData(n, q, p, p, D3)    
    fit <- mqr(mydata$Y, mydata$X)
    D3hat <- fit$Dnew
	D2hat <- TransferModalUnfoldings(D3hat,3,2,p,K,q)
    
    # Example 2
    # The usage of function "mqr_dr()"	
    fit_dr <- mqr_dr(mydata$Y, mydata$X)
    D3hat <- fit_dr$Dnew
	D2hat <- TransferModalUnfoldings(D3hat,3,2,p,K,q)	
    opt <- fit_dr$rk_opt
 
    # Example 3 
    # The usage of function "mvrblockwise()"
    n <- 200
    q <- 5
    s <- 3
    p <- 100
    B <- matrix(runif(q*s, 2,3), s)
    X <- matrix(rnorm(n*p),n,p)
    Y <- X[,1:s]%*%B + matrix(rnorm(n*q),n)
    fit <- mvrblockwise(Y,X) #See details in the function "mvrblockwise"
    fit$activeX
    fit$Bhat
    which(rowSums(fit$Bhat^2)>0)
    fit$muhat
    
    # Example 4
    # The usage of function "mvrcolwise()"
    n <- 200
    q <- 5
    s <- 3
    p <- 100
    B <- matrix(runif(q*s, 2,3), s)
    X <- matrix(rnorm(n*p),n,p)
    Y <- X[,1:s]%*%B + matrix(rnorm(n*q),n)
    fit <- mvrcolwise(Y,X) #See details in the function "mvrcolwise"
    fit$activeX
    fit$Bhat
    which(rowSums(fit$Bhat^2)>0)
    fit$muhat 
	
# References

Symmetric Tensor Estimation for Quadratic Regression. Manuscript.

# Development
The R-package is developed by Xu Liu (liu.xu@sufe.edu.cn).
