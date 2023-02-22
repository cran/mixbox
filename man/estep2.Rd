\name{estep2}
\alias{estep2}
\title{Computing posteriors expected value.}
\description{For a finite mixture model with density function
\deqn{
{\cal{M}}(\bold{y}|\bold{\Psi})=\sum_{g=1}^{G} \omega_{g} f_{\bold{Y}}(\bold{y}, \bold{\Theta}_g),
}
we use method of Basford et al. (1997) for Computing observed Fisher information matrix. Based on above representation for density function of a finite mixture model, we can write
\deqn{
\bold{Y}| {T}, W \sim N_{d}\bigl(\bold{\mu}+\bold{\lambda}{T}, {W} \Sigma\bigr), {T}| W \sim HN({0}, {W}), W  \sim f_W(w| \bold{\theta}).
}
The required posteriors expectations are \eqn{E\bigl(W^{-1}|\bold{y},\bold{\Psi}\bigr)}, \eqn{E\bigl(W^{-1}T|\bold{y},\bold{\Psi}\bigr)}, and \eqn{E\bigl(W^{-1}T^2|\bold{y},\bold{\Psi}\bigr)}.
}
\usage{estep2(Y, G, weight, mu, sigma, lambda, family, skewness, param, theta, tick, h, N, PDF)}
\arguments{
\item{Y}{an \eqn{n\times d} matrix of observations. }
\item{G}{the number of components. }
\item{weight}{a vector of weight parameters (or mixing proportions). }
\item{mu}{a list of location vectors of \code{G} components. }
\item{sigma}{a list of dispersion matrices of \code{G} components. }
\item{lambda}{a list of skewness matrices of \code{G} components. }
\item{family}{name of the mixing distribution. By default \code{family = "constant"} that corresponds to the finite mixture of multivariate normal (or skew normal) distribution. Other candidates for family name are: "bs" (for Birnbaum-Saunders), "burriii" (for Burr type iii), "chisq" (for chi-square), "exp" (for exponential), "f" (for Fisher), "gamma" (for gamma), "gig" (for generalized inverse-Gaussian), "igamma" (for inverse-gamma), "igaussian" (for inverse-Gaussian), "lindley" (for Lindley), "loglog" (for log-logistic), "lognorm" (for log-normal), "lomax" (for Lomax), "ptstable" (for positive \eqn{\alpha}-stable), "ptstable" (for polynomially tilted \eqn{\alpha}-stable), "rayleigh" (for Rayleigh), and "weibull" (for Weibull). }
\item{skewness}{logical statement.If \code{skewness = "FALSE"} means that a symmetric model is fitted to each component (cluster). If \code{skewness = "TRUE"} a skewed model is fitted to each component. }
\item{param}{name of the elements of \eqn{{\bold{\theta)}}} as the parameter vector of mixing distribution with density function \eqn{f_W(w|{\bold{\theta)}}}. }
\item{theta}{a list of maximum likelihood estimator for \eqn{{\bold{\theta)}}} across \code{G} components. }
\item{tick}{ a binary vector whose length depends on type of family. The elements of \code{tick} are either \code{0} or \code{1}. If element of \code{tick} is \code{0}, then the correspoding element of \eqn{\bold{\theta}} is not considered in the formula of \eqn{f_W(w|{\bold{\theta)}}} for computing the required posterior expectations. If element of \code{tick} is \code{1}, then the corresponding element of \eqn{\bold{\theta}} is considered in the formula of \eqn{f_W(w|{\bold{\theta)}}}. For instance, if \code{family = "gamma"} and either its shape or rate parameter is one, then \code{tick = c(1)}. This is while, if \code{family = "gamma"} and both of the shape and rate parameters are in the formula of \eqn{f_W(w|{\bold{\theta)}}}, then \code{tick = c(1, 1)}.}
\item{h}{a positive small value for computing numerical derivative of \eqn{f_W(w| {\bold{\theta)}})} with respect to \eqn{{\bold{\theta)}}}, that is
\eqn{\partial/ \partial {\bold{\theta)}} f_W(w| {\bold{\theta)}}}. }
\item{N}{a large integer number for computing Monte Carlo approximation of expected value within the E-step of the EM algorithm. }
\item{PDF}{expression for mixing density function \eqn{f_W(w| {\bold{\theta)}}}. }
}

\value{
The required posteriors expectations for approximating the observed Fisher information matrix for canonical or unrestricted finite mixture model. These include
\eqn{\tau_{ig}=E\bigl( Z_{ig}=1 \vert {\bold{y}}_{i}, \hat{\bold{\Psi}} \bigr) = {\omega}_g f_{\bold{Y}}\bigl({\bold{y}}_{i}\big|{\bold{\Theta}_g}\bigr)/\bigl[ \sum_{g=1}^{G} {{\omega}}_g f_{\bold{Y}}\bigl({\bold{y}}_{i}\big|{\bold{\Theta}_g}\bigr) \bigr]},
\eqn{E\bigl( W^{-1}\big \vert {\bold{y}}_{i}, \hat{\bold{\Psi}} \bigr)},
\eqn{E\bigl( {\bold{U}}W^{-1}\big \vert {\bold{y}}_{i}, \hat{\bold{\Psi}} \bigr)},
\eqn{E\bigl( {\bold{U}}{\bold{U}}^{\top}W^{-1}\big \vert {\bold{y}}_{i}, \hat{{\bold{\Psi}}} \bigr)}, and
\eqn{E\bigl( \partial f_W(w| {\bold{\theta)}})/\partial {\bold{\theta)}} \big \vert {\bold{y}}_{i}, \hat{{\bold{\Psi}}} \bigr)}
}
\references{
K. E. Basford, D. R. Greenway, G. J. McLachlan, and D. Peel, (1997). Standard errors of fitted means under normal mixture, \emph{Computational Statistics}, 12, 1-17.
}
\author{Mahdi Teimouri}
\examples{
\donttest{
      n <- 100
      G <- 2
 weight <- rep( 0.5, 2 )
    mu1 <- rep(-5  , 2 )
    mu2 <- rep( 5  , 2 )
 sigma1 <- matrix( c( 0.4, -0.20, -0.20, 0.5 ), nrow = 2, ncol = 2 )
 sigma2 <- matrix( c( 0.5,  0.20,  0.20, 0.4 ), nrow = 2, ncol = 2 )
lambda1 <- diag( c(-5, -5) )
lambda2 <- diag( c( 5,  5) )
 theta1 <- c( 10 )
 theta2 <- c( 20 )
     mu <- list( mu1, mu2 )
  sigma <- list( sigma1 , sigma2 )
 lambda <- list( lambda1, lambda2)
  theta <- list( theta1 , theta2 )
  param <- c( "a" )
    PDF <- quote( (a/2)^(a/2)*x^(-a/2 - 1)/gamma(a/2)*exp( -a/(2*x) ) )
  tick  <- rep( 1, 2 )
theta10 <- c( 10, 10 )
theta20 <- c( 20, 20 )
 theta0 <- list( theta10 , theta20 )
      Y <- rmix( n, G, weight, model = "unrestricted", mu, sigma, lambda,
      family = "igamma", theta0)
    out <- estep2( Y[, 1:2], G, weight, mu, sigma, lambda, family = "igamma",
    skewness = "TRUE", param, theta, tick, h = 0.001, N = 3000, PDF)
}
}
\keyword{internal}
