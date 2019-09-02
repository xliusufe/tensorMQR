\name{generateData}
\alias{generateData}
\docType{package}
\title{
  Generate data from MQR model.
}
\description{
  Generate data for a high-dimensional multiresponse quadratic regression, with or without aparsity assumptions.
}

\usage{
generateData(n, q, s, p, D3,SigmaX=diag(p-1),sigma2=0.2,seed_id, t=0.0, rho=0.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Sample size.}
  
  \item{q}{The number of responses, \eqn{q\geq1}.}
  
  \item{s}{The true covariates associating to response, \eqn{s\geq1}.}
    
  \item{p}{The number of covariates, \eqn{p\geq1}.}
  
  \item{D3}{The mode of unfolding \eqn{D_{(3)}}.}
  
  \item{SigmaX}{Covariance of \eqn{X}. Default is identity matrix.}

  \item{sigma2}{err variance. Default is 0.1.}
  
  \item{seed_id}{Seed of generator.}
  
  \item{rho}{The correlation of \eqn{\epsilon_j} and \eqn{\epsilon}_k, where \eqn{j,k\in \{1,\cdots,q\}}}.
}

\details{This function gives \eqn{qp(p+1)/2} coefficients' estimators of MQR. The core
  tensor is a \eqn{r_1\times r_2\times r_3}-tensor. We choose \eqn{r_1}, \eqn{r_2} 
  and \eqn{r_3}  by \code{BIC} or \code{CV}.
}

\value{
  \item{Y}{Response, a \eqn{n\times q}-matrix.}
  
  \item{X}{Design matrix, a \eqn{n\times p}-matrix.}
}

\references{
  Symmetric Tensor Estimation for Quadratic Regression.
}
\keyword{High-dimensional; Gradient descent; Stiefel manifold;  Tensor estimation; Tucker decomposition. }
\seealso{
  mam_sparse
}
\examples{
  # Example 1
  
  D3 <- matrix(runif(72, 0.7, 1), 2, 36)
  mydata <- generateData(200, 2, 6, 6, D3)
    
  Y <- mydata$Y  
  X <- mydata$X
    
    
  # Example 2
  n <- 500
  p <- 10
  q <- 10
  s <- 7
  s0 <- s
  r10=r20=r30=2
  S3 <- matrix(runif(r10*r20*r30,3,7),nrow = r30)
  T1 <- matrix(rnorm(s0*r10),nrow = s0)
  U1 <- qr.Q(qr(T1))
  T1 <- matrix(rnorm(q*r30),nrow = q)
  U3 <- qr.Q(qr(T1))
  D3 <- U3\%*\%S3\%*\%t(kronecker(U1,U1))
  mydata <- generateData(n,q,s0,p,D3)
}
