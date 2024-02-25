\name{rmix}
\alias{rmix}
\title{Generating realization from finite mixture models.}
\description{The density function of a restricted \eqn{G}-component finite mixture model can be represented as
\deqn{
{\cal{M}}(\bold{y}|\bold{\Psi})=\sum_{g=1}^{G} \omega_{g} f_{\bold{Y}}(\bold{y}, \bold{\Theta}_g),
}
where positive constants \eqn{\omega_{1}, \omega_{2},\cdots,\omega_{G}} are called weight (or mixing proportions) parameters with this properties that \eqn{\sum_{g=1}^{G}\omega_{g}=1} and \eqn{\bold{\Psi} = \bigl(\bold{\Theta}_{1},\cdots, \bold{\Theta}_{G}\bigr)^{\top}} with \eqn{\bold{\Theta}_g=\bigl({\bold{\omega}}_g, {\bold{\mu}}_g, {{\Sigma}}_g, {\bold{\lambda}}_g\bigr)^{\top}}. Herein, \eqn{f_{\bold{Y}}(\bold{y}, \bold{\Theta}_g)} accounts for the density function of random vector \eqn{\bold{Y}} within \eqn{g}-th component that admits the representation given by
\deqn{
{\bf{Y}} \mathop=\limits^d {\bold{\mu}}_{g}+\sqrt{W}{\bold{\lambda}}_{g}\vert{Z}_0\vert + \sqrt{W}{\Sigma}_{g}^{\frac{1}{2}} {\bf{Z}}_1,
}
where \eqn{ {\bold{\mu}}_{g} \in {R}^{d} } is location vector, \eqn{ {\bold{\lambda}}_{g} \in {R}^{d} } is skewness vector, and \eqn{\Sigma_{g}} is a positive definite symmetric dispersion matrix for \eqn{g=1,\cdots,G}. Further, \eqn{W} is a positive random variable with mixing density function \eqn{f_W(w| \bold{\theta}_{g})}, \eqn{ {Z}_0\sim N(0, 1) }, and \eqn{ {\bold{Z}}_1\sim N_{d}\bigl( {\bold{0}}, \Sigma_{g}\bigr) }. We note that \eqn{W}, \eqn{Z_0}, and \eqn{ {\bf{Z}}_1} are mutually independent.
}
\usage{rmix(n, G, weight, model = "restricted", mu, sigma, lambda, family = "constant",
        theta = NULL) }
\arguments{
\item{n}{number of realizations. }
\item{G}{number of components. }
\item{weight}{a vector of weight parameters (or mixing proportions). }
\item{model}{It must be \code{"canonical"}, \code{"restricted"}, or \code{"unrestricted"}. By default \code{model="restricted"}. }
\item{mu}{a list of location vectors of \code{G} components. }
\item{sigma}{a list of dispersion matrices of \code{G} components. }
\item{lambda}{a list of skewness vectors of \code{G} components. If mixture model is symmetric, then a vector of zeros of appropriate size should be considered for the skewness vector of the corresponding component.}
\item{family}{name of mixing distribution. By default \code{family = "constant"} that corresponds to the finite mixture of multivariate normal (or skew normal) distribution. Other candidates for family name are: "bs" (for Birnbaum-Saunders), "burriii" (for Burr type iii), "chisq" (for chi-square), "exp" (for exponential), "f" (for Fisher), "gamma" (for gamma), "gigaussian" (for generalized inverse-Gaussian), "igamma" (for inverse-gamma), "igaussian" (for inverse-Gaussian), "lindley" (for Lindley), "loglog" (for log-logistic), "lognorm" (for log-normal), "lomax" (for Lomax), "pstable" (for positive \eqn{\alpha}-stable), "ptstable" (for polynomially tilted \eqn{\alpha}-stable), "rayleigh" (for Rayleigh), and "weibull" (for Weibull). }
\item{theta}{a list of maximum likelihood estimator(s) for \eqn{\theta} (parameter vector of mixing distribution) across \code{G} components. By default it is \code{NULL}. }
}
%\details{}
\value{a matrix with \eqn{n} rows and \eqn{d + 1} columns. The first \eqn{d} columns constitute \eqn{n} realizations from random vector \eqn{\bold{Y}=(Y_1,\cdots,Y_d)^{\top}} and the last column is the label of realization \eqn{\bold{Y}_i} ( for \eqn{ i = 1, \cdots n} ) indicating the component that \eqn{\bold{Y}_i} is coming from.}
\author{Mahdi Teimouri}
\examples{
\donttest{
 weight <- rep( 0.5, 2 )
    mu1 <- rep(-5  , 2 )
    mu2 <- rep( 5  , 2 )
 sigma1 <- matrix( c( 0.4, -0.20, -0.20, 0.4 ), nrow = 2, ncol = 2 )
 sigma2 <- matrix( c( 0.4,  0.10,  0.10, 0.4 ), nrow = 2, ncol = 2 )
lambda1 <- matrix( c( -4, -2,  2,  5 ), nrow = 2, ncol = 2 )
lambda2 <- matrix( c(  4,  2, -2, -5 ), nrow = 2, ncol = 2 )
 theta1 <- c( 10, 10 )
 theta2 <- c( 20, 20 )
     mu <- list( mu1, mu2 )
  sigma <- list( sigma1 , sigma2 )
 lambda <- list( lambda1, lambda2)
  theta <- list( theta1 , theta2 )
      Y <- rmix( n = 100, G = 2, weight, model = "canonical", mu, sigma, lambda,
           family = "igamma", theta )
}
}
