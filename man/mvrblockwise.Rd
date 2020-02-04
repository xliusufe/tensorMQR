\name{mvrblockwise}
\alias{mvrblockwise-function}
\alias{mvrblockwise}
\docType{package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate coefficients of high-dimensional multivariate regression for the grouped-column-wise
}
\description{
This function provides the coefficient matrix estimator of high-dimensional multivariate regression (MVR) with penalty \code{LASSO}, 
\code{MCP} or \code{SCAD}). The tuning parameter is selected by \code{BIC} (the default), \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}.
}
\usage{
mvrblockwise(Y, X ,Z=NULL, method="BIC", ncv=10, penalty="LASSO", isPenColumn=TRUE, 
             group=NULL, lambda=NULL, nlam=50, intercept=TRUE, lam_min=1e-4, 
             eps=1e-6, max_step=50, gamma_pen=2, dfmax=NULL, alpha=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{The response, a vector of size \eqn{n} or  a matrix of size \eqn{n\times q}.}
  
  \item{X}{The covariates to be penalized, a matrix with dimension \eqn{n\times p}. }

  \item{Z}{The covariates without penalization, a matrix with dimension \eqn{n\times d}. The default is \code{NULL}. }
  
  \item{method}{The method to be applied to select parameters.  Either \code{BIC}
                (the default), \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}.}
    
  \item{ncv}{The number of cross-validation folds.  Default is 10. If \code{method} is not \code{CV}, \code{ncv} is useless.}  
  
  \item{penalty}{The penalty to be applied to the model. Either \code{LASSO} (the default), \code{MCP} or \code{SCAD}.}
 
  \item{isPenColumn}{A logical value indicating whether the coefficients associating with \eqn{X_j} that affects whole response 
                    \eqn{Y} is penalized.  Default is \code{TRUE}. If \code{isPenColumn} is \code{TRUE}, the coefficients associating with
                    \eqn{X_j} that affects simultaneously whole response \eqn{y} is penalized for each \eqn{j\in \{1,\cdots,p\}}. 
                    If \code{isPenColumn} is \code{FALSE}, the coefficients associating with \eqn{X_j} that affects single response
                    \eqn{Y_l} is penalized for each \eqn{j\in \{1,\cdots,p\}}, where \eqn{l\in \{1,\cdots,q\}}.} 
   
  \item{group}{A integer vector describing the grouping of the coefficients. For example, we can preset \code{group = rep(1:G,each=K)}. 
               If no grouping, \code{group = rep(1:ncol(X))}. The default is \code{group = rep(1:ncol(X))}.}      
  
  \item{lambda}{A user-specified sequence of lambda values.  By default,
        a sequence of values of length \code{nlam} is computed, equally
        spaced on the log scale.}
  
  \item{nlam}{The number of lambda values. Default is 50.}
  
  \item{intercept}{Should intercept(s) be fitted (default=\code{TRUE}) or set to zero (\code{FALSE})?}  
   
  \item{lam_min}{The smallest value for lambda, as a fraction of
                 lambda.max.  Default is 1e-3.}
                 
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
              relative change in any coefficient is less than \code{eps1}.  
              Default is \code{1e-4}.}
              
  \item{max_step}{Maximum number of iterations. Default is \code{50}.}
                  
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

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
   \item{Bhat}{Estimator of coefficients of \eqn{X}.}
   
   \item{rss}{Residual sum of squares (RSS).}
   
   \item{activeX}{The active set of \eqn{X}. It is a \eqn{p} dimensional vector.}
  
   \item{lambda}{The sequence of regularization parameter values in the path.}
   
   \item{selectedID}{The index of \code{lambda} corresponding to
                    \code{lambda_opt}.}   
                    
   \item{lambda_opt}{The value of \code{lambda} with the minimum
                    \code{BIC} value.}   
                    
   \item{bic}{BIC value used to select variables.}  
   
   \item{muhat}{Estimator of intercept \eqn{\mu}. It is \code{NULL} if \code{intercept} is \code{FALSE}.}
   
   \item{Chat}{Estimator of coefficients of \eqn{Z}. \code{Chat} is \code{NULL} if \code{Z} is \code{NULL}.}
   
   \item{group}{The input group.}
                    
   \item{Y}{Response \eqn{Y}.}
  
   \item{X}{Design matrix \eqn{X}.}                    
}

\references{
Symmetric Tensor Estimation for Quadratic Regression. Manuscript.
}

\examples{
	library(tensorMQR1)
	
	#example 1
	n <- 200
	q <- 5
	s <- 3
	p <- 100
	B <- matrix(runif(q*s, 2,3), s)
	X <- matrix(rnorm(n*p),n,p)
	Y <- X[,1:s]\%*\%B + matrix(rnorm(n*q),n)
	fit <- mvrblockwise(Y,X)
	fit$activeX
	fit$Bhat
	which(rowSums(fit$Bhat^2)>0)
	fit$muhat
	
	#example 2
	K = 5
	n <- 200
	q <- 5
	s <- 4
	p <- 100
	B1 <- matrix(runif(q*K, 2,3), K)
	B2 <- matrix(0,2*K,q)
	B3 <- matrix(runif(q*(s-1)*K, 2,3), (s-1)*K)
	B <- rbind(B1,B2,B3)
	X <- matrix(rnorm(n*p*K),n)
	Y <- X[,1:((s+2)*K)]\%*\%B + matrix(rnorm(n*q),n)
	group <- rep(1:p,each=K)
	fit <- mvrblockwise(Y,X,group=group,isPenColumn=TRUE)
	which(fit$activeX==1)
	fit$Bhat
	which(rowSums(fit$Bhat^2)>0)
	fit$muhat	
	
	#example 3
	K = 5
	n <- 200
	q <- 5
	s <- 4
	d <- 3
	p <- 100
	B1 <- matrix(runif(q*K, 2,3), K)
	B2 <- matrix(0,2*K,q)
	B3 <- matrix(runif(q*(s-1)*K, 2,3), (s-1)*K)
	B <- rbind(B1,B2,B3)
	C <- matrix(runif(q*d, 1,2), d)
	X <- matrix(rnorm(n*p*K),n)
	Z <- matrix(rnorm(n*d),n)
	Y <- X[,1:((s+2)*K)]\%*\%B + Z\%*\%C + matrix(rnorm(n*q),n)
	group <- rep(1:p,each=K)
	fit <- mvrblockwise(Y,X,Z,group=group,isPenColumn=TRUE)
	which(fit$activeX==1)
	fit$Bhat
	which(rowSums(fit$Bhat^2)>0)
	fit$Chat
	fit$muhat		
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Multivariate regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ Variable selection }% __ONLY ONE__ keyword per line
