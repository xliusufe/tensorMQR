\name{TransferModalUnfoldings}
\alias{TransferModalUnfoldings-function}
\alias{TransferModalUnfoldings}
\docType{package}
\title{Transfer a tensor's modal unfoldings to another.}
\description{
  Transfer a tensor's modal unfoldings to another.
}

\usage{
TransferModalUnfoldings(S, d1, d2 , r1, r2, r3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
  A mode-d1-unfolding of a tensor with size \eqn{r_1\times r_2\times r_3}, input unfolding
}
  \item{d1}{
  An integer, the mode of unfolding \eqn{S_{(d_1)}}
}
  \item{d2}{
   An integer, the mode of output unfolding \eqn{S_{(d_2)}}
}
  \item{r1}{
  The fist dimension of tensor 
}
  \item{r2}{
  The second dimension of tensor
}
  \item{r3}{
  The third dimension of tensor 
}
}

\details{
  This function transfers an input mode-d1-unfolding \eqn{S_{(d_1)}} 
  to mode-d2-unfolding \eqn{S_{(d_2)}}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{D}{The output mode-d2-unfolding, \eqn{S_{(d_2)}}}
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
  D1 <- matrix(1:24,nrow = 4) # A tensor unfolding with size 4*6
  D2 <- TransferModalUnfoldings(D1,1,2,4,3,2)
}
