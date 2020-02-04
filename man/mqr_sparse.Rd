\name{mqr_sparse}
\alias{mqr_sparse-function}
\alias{mqr_sparse}
\docType{package}
\title{
  Fit MQR with sparsity assumption and fixed ranks.
}
\description{
  Fit a high-dimensional multiresponse quadratic regression with or without aparsity assumptions, and given ranks given ranks \eqn{r_1, r_2, r_3}. The multivariate sparse group \code{lasso} (\code{mcp} or \code{scad}) and the steepest gradient descent algorithm are used to estimate tensor for sparsity situation. The tuning parameter is selected by \code{BIC} (the default), \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}.
}

\usage{
mqr_sparse(Y, X, r1 = NULL, r3 = NULL, method = "BIC", ncv = 10, isPenU = 0, isPenColumn = 1, 
           penalty = "LASSO", lambda = NULL, SUV = NULL,  nlam = 50,isSym=TRUE, 
           initMethod="LASSO", lam_min = 1e-4, ftol = 1e-6, maxstep = 20,  maxstep1 = 20, 
           eps = 1e-4,  thresh=1e-4, gamma_pen = 2,  dfmax = NULL, alpha = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times q} numeric design matrix for the model.}
  
  \item{r1}{The first dimension of single value matrix of the tensor. Default is 2.}
  
  \item{r3}{The third dimension of single value matrix of the tensor. Default is 2.}
  
  \item{method}{The method to be applied to select parameters.  Either \code{BIC} (the default), \code{AIC}, \code{EBIC}, \code{CV}, 
                or \code{GCV}.}
                
  \item{ncv}{The number of cross-validation folds.  Default is 10. If \code{method} is not \code{CV}, \code{ncv} is useless.}
  
  \item{isPenU}{A logical value indicating whether the rows of \eqn{U} is penalized.  Default is \code{FALSE}. If \code{isPenU} 
                is \code{FALSE}, the coefficients associating with \eqn{X_j} is penalized for each \eqn{j\in \{1,\cdots,p\}}.}
  
  \item{isPenColumn}{A logical value indicating whether the coefficients associating with \eqn{X_j} that affects whole response 
                    \eqn{y} is penalized.  Default is \code{TRUE}. If \code{isPenU} is \code{TRUE}, the coefficients associating with                     
                    \eqn{X_j} that affects whole response \eqn{y} is penalized for each \eqn{j\in \{1,\cdots,p\}}. If \code{isPenU} is                    
                    \code{FALSE}, the coefficients associating with \eqn{X_j} that affects single response \eqn{y_l} is penalized for each                
                    \eqn{j\in \{1,\cdots,p\}}, where \eqn{l\in \{1,\cdots,q\}}.}
  
  \item{penalty}{The penalty to be applied to the model. Either "LASSO" (the default),  
                 "SCAD", or "MCP".}
        
  \item{lambda}{A user-specified sequence of lambda values.  By default,
        a sequence of values of length \code{nlam} is computed, equally
        spaced on the log scale.}
        
  \item{SUV}{A user-specified list of initial coefficient matrix of \eqn{S},
              \eqn{U}, \eqn{V}. By default, initial matrices are provided randomly.}
   
  \item{isSym}{A logical value indicating whether restrict tensor to be symmetric. If isSym is \code{TRUE} (the default), 
             we decompose the tensor to be \eqn{D_{(3)}=VS_{(3)}(U\otimes U)^T}, where the core tensor \eqn{S} is symmetric, 
             and both \eqn{U} and \eqn{V} belong to Stiefel manifold. If isSym is \code{FALSE}, we decompose 
             the tensor to \eqn{S}, \eqn{A}, \eqn{B}, \eqn{C}, that is \eqn{D_{(3)}=CS_{(3)}(B\otimes A)^T}, 
             where the tensor is treated  as being asymmetric.}
  
  \item{initMethod}{One can estimate the initial tensor \eqn{D_{(3)}} as a metrix by choosing a penalty to penalize group-column wise. \code{initMethod}
                    can be \code{LASSO}, \code{MCP} or \code{SCAD}. The default is \code{LASSO}.}     
   
  \item{nlam}{The number of lambda values. Default is 20.}
  
  \item{lam_min}{The smallest value for lambda, as a fraction of
                 lambda.max.  Default is 1e-3.}
                 
  \item{ftol}{Convergence threshhold for the Curvilinear search.  The algorithm iterates until the
              relative change in any coefficient is less than \code{eps1}.  
              Default is \code{1e-6}.}
              
  \item{maxstep}{The maximum iterates number of the steepest gradient descent method.  
                  Default is \code{20}.} 
                  
  \item{max_step1}{Maximum number of outer iterations.  Default is 20.}
                 
  \item{thresh}{The threshold to numerically determine which coefficients are zeros. Since the steepest projected 
                gradient descent method with the approximated penalty can not shrink the estimated row of true zero row 
                of U to exactly zero, we need to determine a numerical threshold. Default is \code{1e-6}.}
   
  \item{eps}{Convergence threshhold for the outer loop.  The algorithm iterates until the
              relative change in any coefficient is less than \code{eps1}.  
              Default is \code{1e-4}.}
              
  \item{gamma_pen}{The tuning parameter of the MCP/SCAD penalty (see details).}
  
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
  This function gives \eqn{qp(p+1)/2} coefficients' estimators of MQR. The core 
  tensor is a \eqn{r_1\times r_1\times r_3}-tensor. \eqn{r_1} and \eqn{r_3}  are fixed.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{betapath}{Solution path of \eqn{\beta}.}
  
  \item{rss}{Residual sum of squares (RSS).}
  
  \item{df}{Degrees of freedom.}
  
  \item{lambda}{The sequence of regularization parameter values in the path.}
  
  \item{lambda_opt}{The value of \code{lambda} with the minimum
                    \code{BIC} value.}
  
  \item{selectedID}{The index of \code{lambda} corresponding to
                    \code{lambda_opt}.}
  
  \item{activeF}{The active set of \eqn{U}. If \code{isPenColumn} is TRUE, \code{activeF} is same as \code{activeX}}
  
  \item{activeX}{The active set of coefficients associtating with \eqn{X}. If \code{isPenColumn} is TRUE, \code{activeX} is same as \code{activeF}}
  
  \item{Snew}{Estimator of \eqn{S_3}.}
  
  \item{Unew}{Estimator of \eqn{U}.}
  
  \item{Vnew}{Estimator of \eqn{V}.}
  
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
   
  fit <- mqr_sparse(mydata$Y, mydata$X)
  activeX <- fit$activeX
}