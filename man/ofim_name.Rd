\name{ofim_name}
\alias{ofim_name}
\title{Creating name for columns and rows of OFI matrix.}
\description{This function creates name for columns and rows of OFI matrix. }
\usage{ofim_name(G, weight, Dim, lambda, model = "restricted", family = "constant",
skewness = "FALSE", param )}
\arguments{
\item{G}{number of components. }
\item{weight}{a vector of weight parameters (or mixing proportions). }
\item{Dim}{dimension size. }
\item{lambda}{a list of skewness matrices of \code{G} components. }
\item{model}{it must be \code{"canonical"}, \code{"restricted"}, or \code{"unrestricted"}. By default \code{model = "restricted"}.}
\item{family}{name of mixing distribution. By default \code{family = "constant"} that corresponds to the finite mixture of multivariate normal (or skew normal) distribution. Other candidates for family name are: "bs" (for Birnbaum-Saunders), "burriii" (for Burr type iii), "chisq" (for chi-square), "exp" (for exponential), "f" (for Fisher), "gamma" (for gamma), "gigaussian" (for generalized inverse-Gaussian), "igamma" (for inverse-gamma), "igaussian" (for inverse-Gaussian), "lindley" (for Lindley), "loglog" (for log-logistic), "lognorm" (for log-normal), "lomax" (for Lomax), "pstable" (for positive \eqn{\alpha}-stable), "ptstable" (for polynomially tilted \eqn{\alpha}-stable), "rayleigh" (for Rayleigh), and "weibull" (for Weibull). }
\item{skewness}{logical statement. By default \code{skewness = "TRUE"} which means that a skewed model is fitted to each component (cluster). If \code{skewness = "FALSE"}, then a symmetric model is fitted to each component. }
\item{param}{name of the elements of \eqn{\bold{\theta}} as the parameter vector of mixing distribution with density function \eqn{f_W(w; {\bold{\theta})}}. }
}
\value{ Vector of names for columns and rows of OFI matrix. }
\author{Mahdi Teimouri}
%\examples{
%\donttest{}
%}
\keyword{internal}
