# tensorMam
 Symmetric Tensor Estimation for Quadratic regression.
 
  For a high-dimensional multiresponse Quadratic regression (MQR) with or without aparsity assumptions, 
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

    library(tensorMQE)

    D3 <- matrix(runif(108, 0.7, 1), 3, 36)
    mydata <- generateData(200, 3, 6, 6, D3)    
    fit <- mqr(mydata$Y, mydata$X)
    coeff <- fit$Dnew
    
    fit_dr <- mqr_dr(mydata$Y, mydata$X)
    opt <- fit_dr$rk_opt
 
# References

Symmetric Tensor Estimation for Quadratic regression. Manuscript.

# Development
The R-package is developed by Xu Liu (liu.xu@sufe.edu.cn).
