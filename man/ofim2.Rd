\name{ofim2}
\alias{ofim2}
\title{Computing observed Fisher information matrix for unrestricted or canonical finite mixture model.}
\description{This function computes the observed Fisher information matrix for a given unrestricted or canonical finite mixture model. For this, we use the method of Basford et al. (1997). The density function of each \eqn{G}-component finite mixture model is given by
\deqn{
g({\bold{y}}|\Psi)=\sum_{g=1}^{G} \omega_{g} f_{\bold{Y}}({\bold{y}}, \Theta_g),
}
where \eqn{\bold{\Psi} = \bigl(\bold{\Theta}_{1},\cdots, \bold{\Theta}_{G}\bigr)^{\top}} with \eqn{\bold{\Theta}_g=\bigl({\bold{\omega}}_g, {\bold{\mu}}_g, {{\Sigma}}_g, {\bold{\lambda}}_g\bigr)^{\top}}. Herein, \eqn{f_{\bold{Y}}(\bold{y}, \bold{\Theta}_g)} accounts for the density function of random vector \eqn{\bold{Y}} within \eqn{g}-th component that admits the representation given by
\deqn{
{\bold{Y}} \mathop=\limits^d {\bold{\mu}}_g+\sqrt{W}{\bold{\lambda}}_g\vert{Z}_0\vert + \sqrt{W}{\Sigma}_g^{\frac{1}{2}} {\bold{Z}}_1,
}
where \eqn{ {\bold{\mu}}_g \in {R}^{d} } is location vector, \eqn{ {\bold{\lambda}}_g \in {R}^{d} } is skewness vector, \eqn{\Sigma_g} is a positive definite symmetric dispersion matrix for \eqn{g =1,\cdots,G}. Further, \eqn{W} is a positive random variable with mixing density function \eqn{f_W(w| \bold{\theta}_g)}, \eqn{ {Z}_0\sim N(0, 1) }, and \eqn{ {\bold{Z}}_1\sim N_{d}\bigl( {\bold{0}}, \Sigma\bigr) }. We note that \eqn{W}, \eqn{Z_0}, and \eqn{ {\bold{Z}}_1} are mutually independent. For approximating the observed Fisher information matrix of the finite mixture models, we use the method of
Basford et al. (1997). Based on this method, using observations \eqn{{\bold{y}}=({{\bold{y}}_{1},{\bold{y}}_{2},\cdots,{\bold{y}}_{n}})^{\top}}, an approximation of the expected information
\deqn{
-E\Bigl\{ \frac{\partial^2 \log L({\bold{\Psi}}) }{\partial \bold{\Psi} \partial \Psi^{\top} } \Bigr\},
}
is give by the observed information as
\deqn{
\sum_{i =1}^{n} \hat{{\bold{h}}}_{i} \hat{{\bold{h}}}^{\top}_{i},
}
where
\deqn{\hat{\bold{h}}_{i} = \frac{\partial}{\partial \bold{\Psi} } \log L_{i}(\hat{\bold{\Psi} }) }
and
\eqn{ \log L(\hat{\bold{\Psi} })= \sum_{i =1}^{n} \log L_{i}(\hat{\bold{\Psi} })= \sum_{i =1}^{n} \log \Bigl\{ \sum_{g=1}^{G} \widehat{{\omega}}_g f_{\bold{Y}}\bigl({\bold{y}}_{i}|\widehat{\bold{\Theta}_g}\bigr)\Bigr\}.}
Herein \eqn{\widehat{\omega}_g} and \eqn{\widehat{\bold{\Theta}_g}} denote the maximum likelihood estimator of \eqn{\omega_g} and \eqn{\bold{\Theta}_g}, for \eqn{g=1,\cdots,G}, respectively. }

\usage{ ofim2(Y, G, weight, model, mu, sigma, lambda, family = "constant", skewness = "FALSE",
      param = NULL, theta = NULL, tick = NULL, h = 0.001, N = 3000, level = 0.05,
    PDF = NULL )
}
\arguments{
\item{Y}{an \eqn{n \times d} matrix of observations. }
\item{G}{number of components. }
\item{weight}{a vector of weight parameters (or mixing proportions). }
\item{model}{It must be \code{"canonical"} or \code{"unrestricted"}. }
\item{mu}{a list of location vectors of \code{G} components. }
\item{sigma}{a list of dispersion matrices of \code{G} components. }
\item{lambda}{a list of skewness vectors of \code{G} components. If model is either \code{"canonical"} or \code{"unrestricted"}, then skewness vactor must be given in matrix form of appropriate size.}
\item{family}{name of the mixing distribution. By default \code{family = "constant"} that corresponds to the finite mixture of multivariate normal (or skew normal) distribution. Other candidates for family name are: "bs" (for Birnbaum-Saunders), "burriii" (for Burr type iii), "chisq" (for chi-square), "exp" (for exponential), "f" (for Fisher), "gamma" (for gamma), "gig" (for generalized inverse-Gaussian), "igamma" (for inverse-gamma), "igaussian" (for inverse-Gaussian), "lindley" (for Lindley), "loglog" (for log-logistic), "lognorm" (for log-normal), "lomax" (for Lomax), "pstable" (for positive \eqn{\alpha}-stable), "ptstable" (for polynomially tilted \eqn{\alpha}-stable), "rayleigh" (for Rayleigh), and "weibull" (for Weibull). }
\item{skewness}{logical statement. By default \code{skewness = "FALSE"} which means that a symmetric model is fitted to each component (cluster). If \code{skewness = "TRUE"}, then a skewed model is fitted to each component. }
\item{param}{name of the elements of \eqn{\theta} as the parameter vector of mixing distribution with density function \eqn{f_W(w|{\bold{\theta)}}}. By default it is \code{NULL}. }
\item{theta}{a list of maximum likelihood estimator for \eqn{\theta} across \code{G} components. By default it is \code{NULL}. }
\item{tick}{a binary vector whose length depends on type of family. The elements of \code{tick} are either \code{0} or \code{1}. If element of \code{tick} is \code{0}, then the corresponding element of \eqn{\bold{\theta}} is not considered in the formula of \eqn{f_W(w|{\bold{\theta)}}} for computing the required posterior expectations. If element of \code{tick} is \code{1}, then the corresponding element of \eqn{\bold{\theta}} is considered in the formula of \eqn{f_W(w|{\bold{\theta)}}}. For instance, if \code{family = "gamma"} and either its shape or rate parameter is one, then \code{tick = c(1)}. This is while, if \code{family = "gamma"} and both of the shape and rate parameters are in the formula of \eqn{f_W(w|{\bold{\theta)}}}, then \code{tick = c(1, 1)}. By default \code{tick = NULL}. }
\item{h}{a positive small value for computing numerical derivative of \eqn{f_W(w| \bold{\theta})} with respect to \eqn{\bold{\theta}}, that is \eqn{\partial/ \partial \theta f_W(w| \bold{\theta})}. By default \eqn{h = 0.001}. }
\item{N}{an integer number for approximating the posterior expected values within the E-step of the EM algorithm through the Monte Carlo method. By default \eqn{N = 3000}. }
\item{level}{significance level \eqn{\alpha} for constructing \eqn{100(1-\alpha)\%} confidence interval. By default \eqn{\alpha = 0.05}. }
\item{PDF}{mathematical expression for mixing density function \eqn{f_W(w| \bold{\theta})}. By default it is \code{NULL}. }
}
%\details{}
\value{
A two-part list whose first part is the observed Fisher information matrix for finite mixture model.
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
 sigma1 <- matrix( c(0.4, -0.20, -0.20, 0.5 ), nrow = 2, ncol = 2 )
 sigma2 <- matrix( c(0.5,  0.20,  0.20, 0.4 ), nrow = 2, ncol = 2 )
lambda1 <- diag( c( 5, -5 ) )
lambda2 <- diag( c(-5,  5 ) )
     mu <- list( mu1, mu2 )
 lambda <- list( lambda1, lambda2 )
  sigma <- list( sigma1 , sigma2  )
    PDF <- quote( (b/2)^(a/2)*x^(-a/2 - 1)/gamma(a/2)*exp( -b/(x*2) ) )
  param <- c( "a","b")
 theta1 <- c( 10, 12 )
 theta2 <- c( 10, 20 )
  theta <- list( theta1, theta2 )
  tick  <- c( 1, 1 )
         Y <- rmix(n, G, weight, model = "unrestricted", mu, sigma, lambda, family = "igamma",
            theta)
    out <- ofim2(Y[, 1:2], G, weight, model = "unrestricted", mu, sigma, lambda,
            family = "igamma", skewness = "TRUE", param, theta, tick, h = 0.001, N = 3000,
           level = 0.05, PDF)
}
}
