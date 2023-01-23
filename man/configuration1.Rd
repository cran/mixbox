 \name{configuration1}
\alias{configuration1}
\title{Output configuration.}
\description{Output design for standard error of the standard error for parameters of the finite mixture models. }
\usage{configuration1(Y, G, weight, mu, sigma, lambda, family, skewness, param, theta,
ofim1, level)}
\arguments{
\item{Y}{an \eqn{n \times d} matrix of observations. }
\item{G}{the number of components. }
\item{weight}{a vector of weight parameters (or mixing proportions). }
\item{mu}{a list of location vectors of \code{G} components. }
\item{sigma}{a list of dispersion matrices of \code{G} components. }
\item{lambda}{a list of skewness vectors of \code{G} components. }
\item{family}{name of the mixing distribution. By default \code{family = "constant"} that corresponds to the finite mixture of multivariate normal (or skew normal) distribution. Other candidates for family name are: "bs" (for Birnbaum-Saunders), "burriii" (for Burr type iii), "chisq" (for chi-square), "exp" (for exponential), "f" (for Fisher), "gamma" (for gamma), "gig" (for generalized inverse-Gausssian), "igamma" (for inverse-gamma), "igaussian" (for inverse-Gausssian), "lindley" (for Lindley), "loglog" (for log-logistic), "lognorm" (for log normal), "lomax" (for Lomax), "pstable" (for positive \eqn{\alpha}-stable), "rayleigh" (for Rayleigh), "ptstable" (for polynomially tilted \eqn{\alpha}-stable), and "weibull" (for Weibull). }
\item{skewness}{logical statement. By default \code{skewness = "TRUE"} which means that a skewed model is fitted to each component (cluster). If \code{skewness = "FALSE"}, then a symmetric model is fitted to each component. }
\item{param}{name of the elements of \eqn{\bold{\theta}} as the parameter vector of mixing distribution with density function \eqn{f_W(w; \bold{\theta})}. }
\item{theta}{a list of maximum likelihood estimator for \eqn{\bold{\theta}} across \code{G} components. }
\item{ofim1}{the observed Fisher information matrix corresponds to the restricted model. }
\item{level}{significance level \eqn{\alpha} for constructing \eqn{100(1-\alpha)\%} confidence interval. By default \eqn{\alpha = 0.05}. }
}
\value{ designated form for output of parameters and their standard errors. }
%\references{}
\author{Mahdi Teimouri}
\examples{
\donttest{
      n <- 2000
      G <- 2
 weight <- rep( 0.5, 2 )
    mu1 <- rep(-5  , 2 )
    mu2 <- rep( 5  , 2 )
 sigma1 <- matrix( c( 0.4, -0.20, -0.20, 0.5 ), nrow = 2, ncol = 2 )
 sigma2 <- matrix( c( 0.5,  0.20,  0.20, 0.4 ), nrow = 2, ncol = 2 )
lambda1 <- c(-5, -5 )
lambda2 <- c( 5,  5 )
 theta1 <- c( 10, 12 )
 theta2 <- c( 10, 20 )
     mu <- list(mu1, mu2)
  sigma <- list( sigma1 , sigma2 )
 lambda <- list( lambda1, lambda2)
  theta <- list( theta1 , theta2 )
  param <- c("a","b")
    PDF <- quote( (b/2)^(a/2)*x^(-a/2 - 1)/gamma(a/2)*exp( -b/(x*2) ) )
  tick  <- c(1, 1)
      Y <- rmix(n, G, weight, model = "restricted", mu, sigma, lambda, family = "igamma", theta)
  ofim  <- ofim1(Y[, 1:2], G, weight, mu, sigma, lambda, family = "igamma",
  skewness = "TRUE", param, theta, tick, h = 0.001, N = 3000, level = 0.05, PDF)
configuration1(Y[, 1:2], G, weight = weight, mu, sigma, lambda, family = "igamma",
skewness = "TRUE", param, theta, ofim1 = ofim, level = 0.05)
}
}
\keyword{internal}
