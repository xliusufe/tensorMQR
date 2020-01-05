\name{mqr_dr}
\alias{mqr_dr-function}
\alias{mqr_dr}
\docType{package}
\title{
  Fit MQR without sparsity assumption, and with ranks selected by \code{BIC} or \code{CV}.
}
\description{
  Fit a low-dimensional multiresponse quadratic regression without aparsity assumptions and with ranks \eqn{r_1, r_3} selected by \code{BIC} or \code{CV}.
} 

\usage{mqr_dr(Y, X, r1_index = NULL, r3_index = NULL, method = "BIC", ncv = 10, 
               SUV = NULL, isSym = TRUE, eps = 1e-6, max_step = 20, max_step1 = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times p} numeric design matrix for the model.}

  \item{r1_index}{A user-specified sequence of \eqn{r_1} values, where 
                   \eqn{r_1} is the first dimension of single value matrix of the 
                   tensor. Default is 
                   \code{r1_index}\eqn{=1,\cdots,\min(\log(n)\rceil,p)}.
                   }
   
  \item{r3_index}{A user-specified sequence of \eqn{r_3} values, where 
                  \eqn{r_3} is the third dimension of single value matrix of the  
                  tensor. Default is 
                  \code{r3_index}\eqn{=1,\cdots,\min(\log(n)\rceil,q)}.
                  }
                  
  \item{method}{The method to be applied to select parameters.  Either \code{BIC} (the default), \code{AIC}, \code{EBIC}, \code{CV}, 
                or \code{GCV}.}
    
  \item{ncv}{The number of cross-validation folds.  Default is 10. If \code{method} is not 
             \code{"CV"}, \code{ncv} is useless.}
                   
  \item{SUV}{A user-specified list of initial coefficient matrix of \eqn{S},
              \eqn{U}, \eqn{V}, which is a list with values \eqn{S},
              \eqn{U}, \eqn{V}. By default, initial matrices are provided randomly.}

  \item{isSym}{A logical value indicating whether restrict tensor to be symmetric. If isSym is \code{TRUE} (the default),  the core tensor 
             is symmetric, and both \eqn{U} and \eqn{V} belong to Stiefel manifold. If isSym is \code{FALSE}, we decompose 
             tensor to \eqn{S}, \eqn{A}, \eqn{B}, \eqn{C}, where we treat the tensor as an asymmetric tensor.}
                
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
             relative change in any coefficient is less than \code{eps}.  Default
             is \code{1e-6}.}
             
  \item{maxstep}{The maximum iterates number of the steepest gradient descent method.  
                  Default is \code{20}.} 
                  
  \item{max_step1}{Maximum number of outer iterations.  Default is 20.}
}

\details{This function gives \eqn{qp(p+1)/2} coefficients' estimators of MQR. The core tensor is a \eqn{r_1\times r_2\times r_3}-tensor. 
         We choose \eqn{r_1}, \eqn{r_2} and \eqn{r_3}  by \code{BIC} or \code{CV}.
}

\value{
  \item{Dnew}{Estimator of \eqn{D_{(3)}}.}
  
  \item{rss}{Residual sum of squares (RSS).}
  
  \item{rk_opt}{The optimal parametres that slected by \code{BIC} or \code{CV}. It is a vector with length 4, which are selected \eqn{r_1} and \eqn{r_3}.}
  
  \item{selected}{Which \eqn{\lambda} is selection.}
  
  \item{Y}{Response \eqn{Y}.}
  
  \item{X}{Design matrix \eqn{X}.}
}
\references{
Symmetric Tensor Estimation for Quadratic Regression.
}
\keyword{High-dimensional; Gradient descent; Stiefel manifold;  Tensor estimation; Tucker decomposition. }
\seealso{
  mqr, mqr_sparse_dr
}
\examples{ 
  D3 <- matrix(runif(72, 0.7, 1), 2, 36) # tensor with size 6*6*2
  mydata <- generateData(200, 2, 6, 6, D3)
  
  fit <- mqr_dr(mydata$Y, mydata$X)
  D3hat <- fit$Dnew
  opt <- fit$rk_opt
}