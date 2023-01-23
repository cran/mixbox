\name{arrange_sigma}
\alias{arrange_sigma}
\title{Extracting diagonal and upper off-diagonal elements of a squre matrix.}
\description{Suppose \eqn{x} is a \eqn{d \times d} square matrix. This function extracts the diagonal and upper off-diagonal elements of a given square matrix.}
\usage{arrange_sigma(x)}
\arguments{
\item{x}{a squre matrix.}
}
\value{ diagonal and upper off-diagonal elements of \eqn{{x}} that consist of \eqn{{x}_{11},{x}_{12},\cdots,{x}_{1d},{x}_{21},{x}_{22},\cdots,{x}_{2d},\cdots,{x}_{dd} }. }
\author{Mahdi Teimouri}
\examples{
\donttest{
x <- matrix( c( 0.4, -0.20, -0.20, 0.5 ), nrow = 2, ncol = 2 )
arrange_sigma(x)
}
}
\keyword{internal}
