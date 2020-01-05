\name{TransferT2P}
\alias{TransferT2P-function}
\alias{TransferT2P}
\docType{package}
\title{Transfer a tensor's modal unfoldings to coefficients of the multivariate quadritic model.}
\description{
  Transfer a tensor's modal unfoldings to coefficients of the multivariate quadritic model.
}

\usage{
TransferT2P(S, d , p, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
  A mode-\eqn{d}-unfolding of a tensor with size \eqn{p\times pq}, input unfolding
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
  This function transfers an input mode-\eqn{d}-unfolding \eqn{S_{(d)}} 
  to coefficients of the multivariate quadritic model  \code{coef}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{coef}{The coefficients of the multivariate quadritic model. \code{coef} is a vector with length \eqn{qp(p+1)/2}.
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
}
