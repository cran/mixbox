\name{arrange_lambda}
\alias{arrange_lambda}
\title{Extracting diagonal and upper off-diagonal elements of a squre matrix.}
\description{Suppose \eqn{x} is a \eqn{d \times d} square matrix. This function extracts the diagonal and upper off-diagonal elements of a given square matrix.}
\usage{arrange_lambda(x)}
\arguments{
\item{x}{a \eqn{m \times n} matrix.}
}
\value{ a list of addresses for elements of the given matrix \eqn{x} as \eqn{{x}_{11},{x}_{12},\cdots,{x}_{1n},{x}_{21},{x}_{22},\cdots,{x}_{2n},\cdots,{x}_{mn} }. }
\author{Mahdi Teimouri}
\examples{
\donttest{
x <- matrix( c( 0.4, -0.20, -0.20, 0.5 ), nrow = 2, ncol = 2 )
arrange_lambda(x)
}
}
\keyword{internal}
