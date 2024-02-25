\name{dmix}
\alias{dmix}
\title{Approximating the density function of the finite mixture models applied for model-based clustering.}
\description{The density function of a \eqn{G}-component finite mixture model can be represented as
\deqn{
g({\bold{y}}|\Psi)=\sum_{g=1}^{G} \omega_{g} f_{\bold{Y}}({\bold{y}}, \Theta_g),
}
where \eqn{\bold{\Psi} = \bigl(\bold{\Theta}_{1},\cdots, \bold{\Theta}_{G}\bigr)^{\top}} with \eqn{\bold{\Theta}_g=\bigl({\bold{\omega}}_g, {\bold{\mu}}_g, {{\Sigma}}_g, {\bold{\lambda}}_g\bigr)^{\top}}. Herein, \eqn{f_{\bold{Y}}(\bold{y}, \bold{\Theta}_g)} accounts for the density function of random vector \eqn{\bold{Y}} within each component. In the restricted case, \eqn{f_{\bold{Y}}(\bold{y}, \bold{\Theta}_g)} admits the representation given by
\deqn{
{\bold{Y}} \mathop=\limits^d {\bold{\mu}}_{g}+\sqrt{W}{\bold{\lambda}}_{g}\vert{Z}_0\vert + \sqrt{W}{\Sigma}_{g}^{\frac{1}{2}} {\bold{Z}}_1,
}
where \eqn{ {\bold{\mu}}_{g} \in {R}^{d} } is location vector, \eqn{ {\bold{\lambda}}_{g} \in {R}^{d} } is skewness vector, \eqn{\Sigma_{g}} is a positive definite symmetric dispersion matrix for \eqn{g=1,\cdots,G}. Further, \eqn{W} is a positive random variable with mixing density function \eqn{f_W(w| \bold{\theta}_{g})}, \eqn{ {Z}_0\sim N(0, 1) }, and \eqn{ {\bold{Z}}_1\sim N_{d}\bigl( {\bold{0}}, \Sigma_{g}\bigr) }. We note that \eqn{W}, \eqn{Z_0}, and \eqn{ {\bold{Z}}_1} are mutually independent. In the canonical or unrestricted case, \eqn{f_{\bold{Y}}(\bold{y}, \bold{\Theta}_g)} admits the representation as
\deqn{
{\bold{Y}} \mathop=\limits^d {\bold{\mu}}_{g}+\sqrt{W}{\bold{\Lambda}}_{g} \vert\bold{Z}_0\vert + \sqrt{W}{\Sigma}_{g}^{\frac{1}{2}} {\bold{Z}}_1,
}
where \eqn{\bold{\Lambda}_{g}} is the skewness matrix and random vector \eqn{\bold{Z}_0} follows a zero-mean normal random vector truncated to the positive hyperplane \eqn{R^{d}} whose independent marginals have variance unity. We note that in the unrestricted case \eqn{\bold{\Lambda}_{g}} is a \eqn{d \times d} diagonal matrix whereas in the canonical case, it is a \eqn{d\times q} matrix and so, random vector \eqn{\bold{Z}_0} follows a zero-mean normal random vector truncated to the positive hyperplane \eqn{R^{q}}.
}
\usage{dmix(Y, G, weight, model = "restricted", mu, sigma, lambda, family = "constant",
skewness = "FALSE", param = NULL, theta = NULL, tick = NULL, N = 3000, log = "FALSE")
 }
\arguments{
\item{Y}{an \eqn{n\times d} matrix of observations. }
\item{G}{number of components. }
\item{weight}{a vector of weight parameters (or mixing proportions). }
\item{model}{it must be \code{"canonical"}, \code{"restricted"}, or \code{"unrestricted"}. By default \code{model = "restricted"}.}
\item{mu}{a list of location vectors of \code{G} components. }
\item{sigma}{a list of dispersion matrices of \code{G} components. }
\item{lambda}{a list of skewness vectors of \code{G} components. If model is either \code{"canonical"} or \code{"unrestricted"}, then skewness vector must be given in matrix form of appropriate size.}
\item{family}{name of mixing distribution. By default \code{family = "constant"} that corresponds to the finite mixture of multivariate normal (or skew normal) distribution. Other candidates for family name are: "bs" (for Birnbaum-Saunders), "burriii" (for Burr type iii), "chisq" (for chi-square), "exp" (for exponential), "f" (for Fisher), "gamma" (for gamma), "gig" (for generalized inverse-Gaussian), "igamma" (for inverse-gamma), "igaussian" (for inverse-Gaussian), "lindley" (for Lindley), "loglog" (for log-logistic), "lognorm" (for log-normal), "lomax" (for Lomax), "pstable" (for positive \eqn{\alpha}-stable), "ptstable" (for polynomially tilted \eqn{\alpha}-stable), "rayleigh" (for Rayleigh), and "weibull" (for Weibull). }
\item{skewness}{a logical statement. By default \code{skewness = "FALSE"} which means that a symmetric model is fitted to each component (cluster). If \code{skewness = "FALSE"}, then a skewed model is fitted to each component. }
\item{param}{name of the elements of \eqn{\bold{\theta}} as the parameter vector of mixing distribution with density function \eqn{f_W(w| \bold{\theta})}. By default it is \code{NULL}. }
\item{theta}{a list of maximum likelihood estimator for \eqn{\bold{\theta}} (parameter vector of the mixing distribution with density function \eqn{f_W(w| \bold{\theta})}), across \code{G} components. By default it is \code{NULL}. }
\item{tick}{ a binary vector whose length depends on type of family. The elements of \code{tick} are either \code{0} or \code{1}. If element of \code{tick} is \code{0}, then the corresponding element of \eqn{\bold{\theta}} is not considered in the formula of \eqn{f_W(w|{\bold{\theta)}}} for computing the required posterior expectations. If element of \code{tick} is \code{1}, then the corresponding element of \eqn{\bold{\theta}} is considered in the formula of \eqn{f_W(w|{\bold{\theta)}}}. For instance, if \code{family = "gamma"} and either its shape or rate parameter is one, then \code{tick = c(1)}. This is while, if \code{family = "gamma"} and both of the shape and rate parameters are in the formula of \eqn{f_W(w|{\bold{\theta)}}}, then \code{tick = c(1, 1)}. By default \code{tick = NULL}.}
\item{N}{an integer number for approximating the \eqn{ g({\bold{y}}|\Psi) }. By default \eqn{N = 3000}. }
\item{log}{if \code{log = "TRUE"}, then it returns the log of the density function. By default it is \code{log = "FALSE"}.}
}
%\details{}
\value{Monte Carlo approximated values of mixture model density function.}
\author{Mahdi Teimouri}
\examples{
\donttest{
      Y <- c(1, 2)
      G <- 2
 weight <- rep( 0.5, 2 )
    mu1 <- rep(  -5, 2 )
    mu2 <- rep(   5, 2 )
 sigma1 <- matrix( c( 0.4, -0.20, -0.20, 0.5 ), nrow = 2, ncol = 2 )
 sigma2 <- matrix( c( 0.5,  0.20,  0.20, 0.4 ), nrow = 2, ncol = 2 )
lambda1 <- c( 5, -5 )
lambda2 <- c(-5,  5 )
     mu <- list( mu1, mu2 )
  sigma <- list( sigma1 , sigma2 )
 lambda <- list( lambda1, lambda2)
    out <- dmix(Y, G, weight, model = "restricted", mu, sigma, lambda, family =
           "constant", skewness = "TRUE", param = NULL, theta = NULL, tick =
           NULL, N = 3000)
}
}
