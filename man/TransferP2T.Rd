\name{TransferP2T}
\alias{TransferP2T-function}
\alias{TransferP2T}
\docType{package}
\title{Transfer coefficients of the multivariate quadritic model to a tensor's modal unfoldings.}
\description{
  Transfer coefficients of the multivariate quadritic model to a tensor's modal unfoldings.
}

\usage{
TransferP2T(coef, d , p, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coef}{
  The coefficients of the multivariate quadritic model, which is a vector with length \eqn{qp(p+1)/2}
}
  \item{d}{
  An integer, the mode of unfolding \eqn{S_{(d)}}
}

  \item{p}{
  The fist dimension of tensor 
}
  \item{q}{
  The third dimension of tensor
}
}

\details{
  This function transfers coefficients of the multivariate quadritic model  \code{coef} to a mode-\eqn{d}-unfolding \eqn{D_{(d)}} 
  of a tensor.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Dd}{A mode-\eqn{d}-unfolding of a tensor with size \eqn{p\times pq}, input unfolding.
}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}
\references{
Symmetric Tensor Estimation for Quadratic Regression.
}
\keyword{High-dimensional; Gradient descent; Stiefel manifold;  Tensor estimation; Tucker decomposition. }
\examples{
  p <- 4
  q <- 3
  D1 <- NULL
  for(j in 1:q){
    D0 <- matrix(runif(p^2,1,3),p)
    D1 <- cbind(D1,(D0+t(D0))/2)
  }   
  coef <- TransferT2P(D1, 1 , p, q)
  D1 <- TransferP2T(coef, 1 , p, q)
  D2 <- TransferP2T(coef, 2 , p, q)
  coef2 <- TransferT2P(D2, 2 , p, q)
}
