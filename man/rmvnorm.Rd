\name{rmvnorm}
\alias{rmvnorm}
\title{Generating from multivariate normal distribution with location vector \eqn{\bold{\mu}} and covariance matrix \eqn{\Sigma_{d \times d}}.}
\description{Using the well-recognized Cholesky decomposition, this function simulates from the density function of a \eqn{d}-dimensional random vector \eqn{\bold{Y}=(Y_1,\cdots,Y_d)^{T}} following a normal distribution with mean vector \eqn{\bold{\mu}} and covariance matrix \eqn{\Sigma_{d \times d}} is
\deqn{
f_{\bold{Y}}(\bold{y})=\frac{1}{(2 \pi)^{\frac{d}{2}}\vert \Sigma\vert ^{-\frac{1}{2} } } \exp\biggl\{-\frac{(\bold{y}-\bold{\mu})^t\Sigma^{-1}(\bold{y}-\bold{\mu})}{2}\bigg\},
}
}
\usage{rmvnorm(n, Mu, Sigma) }
\arguments{
\item{n}{number of realizations. }
\item{Mu}{location vector. }
\item{Sigma}{covariance (dispersion) matrix. }
}
%\details{}
\value{an \eqn{n \times d} matrix of realizations from multivariate normal distribution with mean vector \eqn{\bold{\mu}} and covariance matrix \eqn{\Sigma_{d \times d}}.}
\author{Mahdi Teimouri}
\examples{
\donttest{
      n <- 100
    Mu  <- rep(0, 2)
 Sigma  <- matrix( c( 2, 0.50, 0.50, 2 ), nrow = 2, ncol = 2 )
 rmvnorm(n, Mu, Sigma)
}
}
\keyword{internal}

