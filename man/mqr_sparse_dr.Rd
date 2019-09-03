\name{mqr_sparse_dr}
\alias{mqr_sparse_dr-function}
\alias{mqr_sparse_dr}
\docType{package}
\title{
  Fit MQR with sparsity assumption and ranks selected by \code{BIC} or \code{CV}.
}
\description{
  Fit a high-dimensional multiresponse quadratic regression with or with aparsity assumptions and ranks selected by \code{BIC} or \code{CV}. The multivariate sparse group \code{lasso} (\code{mcp} or \code{scad}) and the steepest gradient descent algorithm are used to estimate
  tensor for sparsity situation. The tuning parameter is selected by \code{BIC} or \code{CV}, which matchs the method of rank selection.
}

\usage{
mqr_sparse_dr(Y, X, r1_index = NULL, r3_index = NULL, method = "BIC", ncv = 10, 
              penalty = "LASSO", isPenU=0,isPenColumn=1, lambda = NULL, SUV = NULL, 
              nlam = 50, lam_min = 0.001, eps = 1e-6, maxstep = 20, maxstep1 = 20, 
              gamma = 2, dfmax = NULL, alpha = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times q} numeric design matrix for the model.}
  
  \item{r1_index}{A user-specified sequence of \eqn{r_1} values, where 
                  \eqn{r_1} is the first dimension of single value matrix of the  
                  tensor. Default is 
                  \code{r1_index}\eqn{=1,\cdots,\min(\lceil\log(n)\rceil,p)}.}
                  
  \item{r3_index}{A user-specified sequence of \eqn{r_3} values, where 
                  \eqn{r_3} is the third dimension of single value matrix of the  
                  tensor. Default is 
                  \code{r3_index}\eqn{=1,\cdots,\min(\lceil\log(n)\rceil,q)}.}
                  
  \item{method}{The method to be applied to select parameters.  Either "BIC"
                (the default), or "CV".}
    
  \item{ncv}{The number of cross-validation folds.  Default is 10. If \code{method} is \code{BIC}, \code{ncv} is useless.}
  
  \item{penalty}{The penalty to be applied to the model. Either "LASSO" (the default), 
                 "SCAD", or "MCP".}
 
  \item{isPenU}{A logical value indicating whether the rows of \eqn{U} is penalized.  Default is \code{FALSE}. If \code{isPenU} 
                is \code{FALSE}, the coefficients associating with \eqn{X_j} is penalized for each \eqn{j\in \{1,\cdots,p\}}.}
  
  \item{isPenColumn}{A logical value indicating whether the coefficients associating with \eqn{X_j} that affects whole response 
                    \eqn{y} is penalized.  Default is \code{TRUE}. If \code{isPenU} is \code{TRUE}, the coefficients associating with                             \eqn{X_j} that affects whole response \eqn{y} is penalized for each \eqn{j\in \{1,\cdots,p\}}. If \code{isPenU} is                            \code{FALSE}, the coefficients associating with \eqn{X_j} that affects single response \eqn{y_l} is penalized for each                        \eqn{j\in \{1,\cdots,p\}}, where \eqn{l\in \{1,\cdots,q\}}.}
  
  \item{lambda}{A user-specified sequence of lambda values.  By default,
                a sequence of values of length \code{nlam} is computed, equally
                spaced on the log scale.}

  \item{SUV}{A user-specified list of initial coefficient matrix of \eqn{S},
              \eqn{U}, \eqn{V}. By default, initial matrices are provided randomly.}
              
  \item{nlam}{The number of lambda values.  Default is 50.}
  
  \item{lam_min}{The smallest value for lambda, as a fraction of
                 lambda.max.  Default is 1e-2.}
                 
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
              relative change in any coefficient is less than \code{eps1}.  
              Default is \code{1e-6}.}
              
  \item{maxstep}{Maximum number of iterations. Default is 20.}
  
  \item{maxstep1}{The maximum iterates number of the steepest gradient descent method.  
                  Default is \code{20}.}

  \item{thresh}{The threshold to numerically determine which coefficients are zeros. Since the steepest projected 
                gradient descent method with the approximated penalty can not shrink the estimated row of true zero row 
                of U to exactly zero, we need to determine a numerical threshold. Default is \code{1e-6}.}

  \item{gamma}{The tuning parameter of the MCP/SCAD penalty (see details).}
  
  \item{dfmax}{Upper bound for the number of nonzero coefficients.
               Default is no upper bound.  However, for large data sets,
               computational burden may be heavy for models with a large number of
               nonzero coefficients.}
    
  \item{alpha}{Tuning parameter for the Mnet estimator which controls
               the relative contributions from the LASSO, MCP/SCAD penalty and the ridge,
               or L2 penalty.  \code{alpha=1} is equivalent to LASSO, MCP/SCAD penalty,
               while \code{alpha=0} would be equivalent to ridge regression.
               However, \code{alpha=0} is not supported; \code{alpha} may be
               arbitrarily small, but not exactly 0.}

}

\details{
  This function gives \eqn{qp(p+1)/2} coefficients' estimators of MAM. The core 
  tensor is a \eqn{r_1\times r_2\times r_3}-tensor. We choose \eqn{r_1}, \eqn{r_2} 
  and \eqn{r_3}  by \code{BIC} or \code{CV}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  
  \item{rss}{Residual sum of squares (RSS).}
  
  \item{df}{Degrees of freedom.}
  
  \item{activeF}{The active set of \eqn{U}. If \code{isPenColumn} is TRUE, \code{activeF} is same as \code{activeX}}
  
  \item{activeX}{The active set of coefficients associtating with \eqn{X}. If \code{isPenColumn} is TRUE, \code{activeX} 
                 is same as \code{activeF}}
  
  \item{Snew}{Estimator of \eqn{S_3}.}
  
  \item{Unew}{Estimator of \eqn{U}.}
  
  \item{Vnew}{Estimator of \eqn{V}.}
  
  \item{lambda}{The sequence of regularization parameter values in the path.}
  
  \item{selectedID }{The index of \code{lambda} corresponding to \code{lambda_opt}.}
  
  \item{lambda_opt }{The value of \code{lambda} with the minimum \code{BIC} or \code{CV} value.}
  
  \item{RSS}{The values of \code{BIC} or \code{CV}, which is a vector.}
  
  \item{rk_opt}{The optimal parametres that slected by \code{BIC} or \code{CV}. It is a vector with length 4, which are 
                selected \eqn{r_1} and \eqn{r_3}.}
  
  \item{Y}{Response \eqn{Y}.}
  
  \item{X}{Design matrix \eqn{X}.}
}
\references{
Symmetric Tensor Estimation for Quadratic Regression.
}
\keyword{High-dimensional; Gradient descent; Stiefel manifold;  Tensor estimation; Tucker decomposition. }
\seealso{
  mqr_dr, mqr_sparse
}

\examples{
  #Example 1

  D3 <- matrix(runif(72, 0.7, 1), 2, 36) # tensor with size 6*6*2
  mydata <- generateData(200, 2, 6, 6, D3)
  
  fit <- mqr_sparse_dr(mydata$Y, mydata$X)
  S3hat <- fit$Snew
  opt <- fit$rk_opt
}