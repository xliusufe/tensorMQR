\name{mqr}
\alias{mqr-function}
\alias{mqr}
\docType{package}
\title{
  Fit MQR without sparsity assumption and with fixed ranks.
}
\description{
  Fit a low-dimensional multiresponse quadratic regression without aparsity assumptions and  with given \eqn{r_1, r_2, r_3}. 
  The steepest gradient descent algorithm is used to estimate tensor.
}

\usage{mqr(Y, X, r1 = NULL, r3 = NULL, SUV = NULL, isSym = TRUE, 
                         eps = 1e-6, max_step = 20, max_step1=20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times p} numeric design matrix for the model.}
  
  \item{r1}{The first dimension of single value matrix of the tensor. Default is 2.}
  
  \item{r3}{The third dimension of single value matrix of the tensor. Default is 2.}
  
  \item{SUV}{A user-specified list of initial coefficient matrix of \eqn{S},
              \eqn{U}, \eqn{V}. By default, initial matrices are provided randomly.}

  \item{isSym}{A logical value indicating whether restrict tensor to be symmetric. If isSym is \code{TRUE} (the default),  the core tensor 
             is symmetric, and both \eqn{U} and \eqn{V} belong to Stiefel manifold. If isSym is \code{FALSE}, we decompose 
             the tensor to \eqn{S}, \eqn{A}, \eqn{B}, \eqn{C}, where the tensor is treated  to be asymmetric.}
             
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
        loss function change in any coefficient is less than \code{eps}.  
        Default is \code{1e-6}.}
    
  \item{maxstep}{The maximum iterates number of the steepest gradient descent method.  
                  Default is \code{20}.} 
                  
  \item{max_step1}{Maximum number of outer iterations.  Default is 20.}
  
 
}

\details{This function gives \eqn{qp(p+1)/2} coefficients' estimators of MQR. The core tensor is a \eqn{r_1\times r_2\times r_3}-tensor. 
        We fixed \eqn{r_1}, \eqn{r_2} and \eqn{r_3} in the function \code{mqr}, but one can choose \eqn{r_1}, \eqn{r_2} and \eqn{r_3}  
        by \code{BIC} or \code{CV}. See details in function \code{mqr_bic} or \code{mqr_cv}.
}

\value{
  \item{Dnew}{Estimator of \eqn{D_{(3)}}.}
  
  \item{rss }{Residual sum of squares (RSS).}
  
  \item{Y}{Response \eqn{Y}.}
  
  \item{X}{Design matrix \eqn{X}.}
}

\references{
Symmetric Tensor Estimation for Quadratic Regression.
}
\keyword{High-dimensional; Gradient descent; Stiefel manifold;  Tensor estimation; Tucker decomposition. }
\seealso{
  mqr_sparse, mqr_bic, mqr_cv
}
\examples{ 
  D3 <- matrix(runif(72, 0.7, 1), 2, 36) # tensor with size 6*6*2
  mydata <- generateData(200, 2, 6, 6, D3)
    
  fit <- mqr(mydata$Y, mydata$X, r1=4, r3= 2)
  D3hat <- fit$Dnew
}