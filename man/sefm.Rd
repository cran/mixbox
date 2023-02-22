\name{sefm}
\alias{sefm}
\title{Approximating the asymptotic standard error for parameters of the finite mixture models based on the observed Fisher information matrix.}
\description{The density function of each finite mixture model can be represented as
\deqn{
{\cal{M}}(\bold{y}|\bold{\Psi})=\sum_{g=1}^{G} \omega_{g} f_{\bold{Y}}(\bold{y}, \bold{\Theta}_g),
}
where positive constants \eqn{\omega_{1}, \omega_{2},\cdots,\omega_{G}} are called weight (or mixing proportions) parameters with this properties that \eqn{\sum_{g=1}^{G}\omega_{g}=1} and \eqn{\bold{\Psi} = \bigl(\bold{\Theta}_{1},\cdots, \bold{\Theta}_{G}\bigr)^{\top}} with \eqn{\bold{\Theta}_g=\bigl({\bold{\omega}}_g, {\bold{\mu}}_g, {{\Sigma}}_g, {\bold{\lambda}}_g\bigr)^{\top}}. Herein, \eqn{f_{\bold{Y}}(\bold{y}, \bold{\Theta}_g)} accounts for the density function of random vector \eqn{\bold{Y}} within \eqn{g}-th component that admits the representation given by
\deqn{
{\bold{Y}} \mathop=\limits^d {\boldsymbol{\mu}}_g+\sqrt{W}{\boldsymbol{\lambda}}_g\vert{Z}_0\vert + \sqrt{W}{\Sigma}_g^{\frac{1}{2}}{\bold{Z}}_1,
}
where \eqn{{\bold{\mu}_g} \in {R}^{d} } is location vector, \eqn{{\bold{\lambda}_g} \in {R}^{d}} is skewness vector, \eqn{\Sigma_g} is a positive definite symmetric dispersion matrix for \eqn{g=1,\cdots,G}. Further, \eqn{W} is a positive random variable with mixing density function \eqn{f_W(w| \bold{\theta}_g)}, \eqn{{Z}_{0}\sim N({0},1)}, and \eqn{{\bold{Z}}_{1}\sim N_{d}\bigl({\bold{0}}, \Sigma_g\bigr)}. We note that \eqn{W}, \eqn{Z_0}, and \eqn{{\bold{Z}}_{1}} are mutually independent. For approximating the asymptotic standard error for parameters of the finite mixture model based on observed Fisher information matrix, we use the method of
Basford et al. (1997). In fact, the covariance matrix of maximum likelihood (ML) estimator \eqn{\hat{\bold{\Psi}}}, can be approximated by the inverse of the observed information matrix as
\deqn{
\sum_{i = 1}^{n} \hat{{\bold{h}}}_{i} \hat{{\bold{h}}}^{\top}_{i},
}
where
\deqn{ \hat{\bold{h}}_{i} = \frac{\partial}{\partial \bold{\Psi} } \log L_{i}(\hat{\bold{\Psi}}), }
and
\eqn{ \log L(\hat{\bold{\Psi}}) = \sum_{i =1}^{n} \log L_{i}(\hat{\bold{\Psi}}) = \sum_{i =1}^{n} \log \Bigl\{ \sum_{g=1}^{G} \widehat{\omega_{g}} f_{\bold{Y}}\bigl({\bold{y}}_{i}|\widehat{\bold{\Theta}_{g}}\bigr)\Bigr\}}. Herein \eqn{\widehat{\omega}_g} and \eqn{\widehat{\bold{\Theta}_g}}, for \eqn{g=1,\cdots,G}, denote the ML estimator of \eqn{\omega_g} and \eqn{\bold{\Theta}_g}, respectively.
}
\usage{sefm(Y, G, weight, model = "restricted", mu, sigma, lambda, family = "constant",
    skewness = "FALSE", param = NULL, theta = NULL, tick = NULL, h = 0.001, N = 3000,
  level = 0.05, PDF = NULL)}
\arguments{
\item{Y}{an \eqn{n\times d} matrix of observations gives the coordinates of the data points. }
\item{G}{number of components. }
\item{weight}{a vector of weight parameters (or mixing proportions). }
\item{model}{it must be \code{"canonical"}, \code{"restricted"}, or \code{"unrestricted"}. By default \code{model = "restricted".} }
\item{mu}{a list of location vectors of \code{G} components. }
\item{sigma}{a list of dispersion matrices of \code{G} components. }
\item{lambda}{a list of skewness vectors of \code{G} components. If model is either \code{"canonical"} or \code{"unrestricted"}, then skewness vactor must be given in matrix form of appropriate size.}
\item{family}{name of mixing distribution. By default \code{family = "constant"} that corresponds to the finite mixture of multivariate normal (or skew normal) distribution. Other candidates for family name are: "bs" (for Birnbaum-Saunders), "burriii" (for Burr type iii), "chisq" (for chi-square), "exp" (for exponential), "f" (for Fisher), "gamma" (for gamma), "gig" (for generalized inverse-Gaussian), "igamma" (for inverse-gamma), "igaussian" (for inverse-Gaussian), "lindley" (for Lindley), "loglog" (for log-logistic), "lognorm" (for log-normal), "lomax" (for Lomax), "pstable" (for positive \eqn{\alpha}-stable), "ptstable" (for polynomially tilted \eqn{\alpha}-stable), "rayleigh" (for Rayleigh), and "weibull" (for Weibull). }
\item{skewness}{a logical statement. By default \code{skewness = "FALSE"} which means that a symmetric model is fitted to each component (cluster). If \code{skewness = "TRUE"}, then a skewed model is fitted to each component. }
\item{param}{name of the elements of \eqn{\bold{\theta}} as the parameter vector of mixing distribution with density function \eqn{f_W(w| \bold{\theta})}. By default it is \code{NULL}. }
\item{PDF}{mathematical expression for mixing density function \eqn{f_W(w| \bold{\theta})}. By default it is \code{NULL}. }
\item{theta}{a list of maximum likelihood estimator for \eqn{\bold{\theta}} across \code{G} components. By default it is \code{NULL}. }
\item{tick}{ a binary vector whose length depends on type of family. The elements of \code{tick} are either \code{0} or \code{1}. If element of \code{tick} is \code{0}, then the corresponding element of \eqn{\bold{\theta}} is not considered in the formula of \eqn{f_W(w|{\bold{\theta)}}} for computing the required posterior expectations. If element of \code{tick} is \code{1}, then the corresponding element of \eqn{\bold{\theta}} is considered in the formula of \eqn{f_W(w|{\bold{\theta)}}}. For instance, if \code{family = "gamma"} and either its shape or rate parameter is one, then \code{tick = c(1)}. This is while, if \code{family = "gamma"} and both of the shape and rate parameters are in the formula of \eqn{f_W(w|{\bold{\theta)}}}, then \code{tick = c(1, 1)}. By default \code{tick = NULL}.}
\item{h}{a positive small value for computing numerical derivative of \eqn{f_W(w| \bold{\theta})} with respect to \eqn{\bold{\theta}}, that is \eqn{\partial/ \partial \theta f_W(w| \bold{\theta})}. By default \eqn{h = 0.001}. }
\item{N}{an integer number for approximating the posterior expected values within the E-step of the EM algorithm through the Monte Carlo method. By default \eqn{N = 3000}. }
\item{level}{significance level \eqn{\alpha} for constructing \eqn{100(1-\alpha)\%} confidence interval. By default \eqn{\alpha = 0.05}. }
}
\details{Mathematical expressions for density function of mixing distribution \eqn{f_W(w|{\bold{\theta}})}, are given as follows.

{{\code{bs}}} {\deqn{f_W(w|{\bold{\theta}}) = \frac {\sqrt{\frac {w}{\beta}}+\sqrt{\frac {\beta}{w}}}{2\sqrt{2\pi}\alpha w}\exp\Biggl\{-\frac {1}{2\alpha^2} \Bigl[\frac{w}{\beta}+\frac{\beta}{w}-2\Bigr]\Biggr\},} where \eqn{{\bold{\theta}}=(\alpha,\beta)^{\top}}. Herein \eqn{\alpha> 0} and \eqn{\beta> 0} are the first and second parameters of this family, respectively.}

{{\code{burrii}}}{\deqn{f_W(w|{\bold{\theta}}) = \alpha \beta w^{-\beta-1}\bigl( 1+w^{-\beta} \bigr) ^{-\alpha-1},} where \eqn{w>0} and \eqn{{\bold{\theta}}=(\alpha, \beta)^{\top}}. Herein \eqn{\alpha> 0} and \eqn{\beta> 0} are the first and second parameters of this family, respectively.}

{{\code{chisq}}}{\deqn{ f_W(w|{{\theta}}) = \frac{2^{-\frac {\alpha}{2}}}{\Gamma\bigl(\frac{\alpha}{2}\bigr)} w^{\frac{\alpha}{2}-1}\exp\Bigl\{-\frac {w}{2} \Bigr\},}
where \eqn{w>0} and \eqn{{{\theta}}=\alpha}. Herein \eqn{\alpha> 0} is the degrees of freedom parameter of this family.}

{\code{exp}}{\deqn{f_W(w|{{\theta}}) =\alpha \exp \bigl\{-\alpha w\bigr\},} where \eqn{w>0} and \eqn{{{\theta}}=\alpha} where \eqn{\alpha> 0} is the rate parameter of this family.}

{\code{f}}{\deqn{ f_W(w|{\bold{\theta}}) = B^{-1}\Bigl(\frac {\alpha}{2},\frac {\beta}{2}\Bigr)\Bigl( \frac {\alpha}{\beta} \Bigr)^{\frac {\alpha}{2}} w^{\frac {\alpha}{2}-1}\Bigl(1 + \alpha\frac {w}{\beta} \Bigr)^{-\left(\frac {\alpha+\beta}{2} \right)},} where \eqn{w>0} and \eqn{B(.,.)} denotes the ordinary beta function. Herein \eqn{{\bold{\theta}}=(\alpha, \beta)^{\top}} where \eqn{\alpha> 0} and \eqn{\beta> 0} are the first and second degrees of freedom parameters of this family, respectively.}

{\code{gamma}}{\deqn{ f_W(w|{\bold{\theta}}) = \frac {\beta^{\alpha}}{\Gamma(\alpha)} \Bigl( \frac{w}{\beta}\Bigr)^{\alpha-1}\exp\bigl\{ - \beta w \bigr\},} where \eqn{w>0} and \eqn{{\bold{\theta}}=(\alpha, \beta)^{\top}}. Herein \eqn{\alpha> 0} and \eqn{\beta> 0} are the shape and rate parameters of this family, respectively.}

{\code{gigaussian}}{\deqn{ f_W(w|{\bold{\theta}}) =\frac{1}{2{\cal{K}}_{\alpha}(\sqrt{\beta \delta})}\Bigl(\frac{\beta}{\delta}\Bigr)^{\alpha/2}w^{\alpha-1}\exp\biggl\{-\frac{\delta}{2w}-\frac{\beta w}{2}\biggr\},} where \eqn{{\cal{K}}_{\alpha}(.)} denotes the modified Bessel function of the third kind with order index \eqn{\alpha} and \eqn{{\bold{\theta}}=(\alpha, \delta, \beta)^{\top}}. Herein \eqn{-\infty <\alpha <+\infty}, \eqn{\delta> 0}, and \eqn{\beta> 0} are the first, second, and third parameters of this family, respectively.}

{\code{igamma}}{\deqn{ f_W(w|{\bold{\theta}}) = \frac{1}{\Gamma(\alpha)} \Bigl( \frac{w}{\beta}\Bigr)^{-\alpha-1}\exp\Bigl\{ - \frac{\beta}{w} \Bigr\},} where \eqn{w>0} and \eqn{{\bold{\theta}}=(\alpha, \beta)^{\top}}. Herein \eqn{\alpha> 0} and \eqn{\beta> 0} are the shape and scale parameters of this family, respectively.}

{\code{igaussian}}{\deqn{ f_W(w|{\bold{\theta}}) =\sqrt{\frac{\beta}{2 \pi w^3}}\exp\biggl\{-\frac{\beta(w - \alpha)^2}{2\alpha^2 w}\biggr\},} where \eqn{w>0} and \eqn{{\bold{\theta}}=(\alpha, \beta)^{\top}}. Herein \eqn{\alpha>0} is the mean and \eqn{\beta> 0} are the first (mean) and second (shape) parameter of this family, respectively.}

{\code{lidley}}{\deqn{f_W(w|{{\theta}}) =\frac{\alpha^2}{\alpha+1} (1+w)\exp \bigl\{-\alpha w\bigr\},} where \eqn{w>0} and \eqn{{{\theta}}=\alpha} where \eqn{\alpha> 0} is the only parameter of this family.}

{\code{loglog}}{\deqn{ f_W(w|{\bold{\theta}}) =\frac{\alpha}{ \beta^{\alpha}} w^{\alpha-1} \left[ \Bigl( \frac {w}{\beta}\Bigr)^\alpha +1\right] ^{-2},} where \eqn{w>0} and \eqn{{\bold{\theta}}=(\alpha, \beta)^{\top}}. Herein \eqn{\alpha> 0} and \eqn{\beta> 0} are the shape and scale (median) parameters of this family, respectively.}

{\code{lognorm}}{\deqn{ f_W(w|{\bold{\theta}}) = \bigl(\sqrt{2\pi} \sigma w \bigr)^{-1}\exp\biggl\{ -\frac{1}{2}\left( \frac {\log w - \mu}{\sigma}\right) ^2\biggr\},} where \eqn{w>0} and \eqn{{\bold{\theta}}=(\mu, \sigma)^{\top}}. Herein \eqn{-\infty<\mu<+\infty} and \eqn{\sigma> 0} are the first and second parameters of this family, respectively.}

{\code{lomax}}{\deqn{ f_W(w|{\bold{\theta}}) = \alpha \beta \bigl( 1+\beta w\bigr)^{-(\alpha+1)},} where \eqn{w>0} and \eqn{{\bold{\theta}}=(\alpha, \beta)^{\top}}. Herein \eqn{\alpha>0} and \eqn{\beta> 0} are the shape and rate parameters of this family, respectively.}

{\code{rayleigh}}{\deqn{ f_W(w|{{\theta}}) = 2\frac {w}{\beta^2}\exp\biggl\{ -\Bigl( \frac {w}{\beta}\Bigr)^2 \biggr\},} where \eqn{w>0} and \eqn{{{\theta}}=\beta}. Herein \eqn{\beta>0} is the scale parameter of this family.}

{\code{weibull}}{\deqn{ f_W(w|{\bold{\theta}}) = \frac {\alpha}{\beta} \Bigl( \frac {w}{\beta} \Bigr)^{\alpha - 1}\exp\biggl\{ -\Bigl( \frac {w}{\beta}\Bigr)^\alpha \biggr\},} where \eqn{w>0} and \eqn{{\bold{\theta}}=(\alpha, \beta)^{\top}}. Herein \eqn{\alpha>0} and \eqn{\beta> 0} are the shape and scale parameters of this family, respectively.}
The density functions of positive \eqn{\alpha}-stable and polynomially tilted \eqn{\alpha}-stable distribution have no closed form and so are not represented here. In what follows, we give four examples. In the first, second, and third examples, we consider three mixture models including: two-component normal, two-component restricted skew \eqn{t}, and two-component restricted skew sub-Gaussian \eqn{\alpha}-stable (SSG) mixture models are fitted to \code{iris}, \code{AIS}, and \code{bankruptcy} data, respectively. In order to approximate the asymptotic standard error of the model parameters, the ML estimators for parameters of skew \eqn{t} and SSG mixture models have been computed through the \code{R} packages \code{EMMIXcskew} (developed by Lee and McLachlan (2018) for skew \eqn{t}) and \code{mixSSG} (developed by Teimouri (2022) for skew sub-Gaussian \eqn{\alpha}-stable). To avoid running package \code{mixSSG}, we use the ML estimators correspond to \code{bankruptcy} data provided by Teimouri (2022). The package \code{mixSSG} is available at \url{https://CRAN.R-project.org/package=mixSSG}. In the fourth example, we apply a three-component generalized hyperbolic mixture model to \code{Wheat} data. The ML estimators of this mixture model have been obtained using the \code{R} package \code{MixGHD} available at \url{https://cran.r-project.org/package=MixGHD}. Finally, we note that if parameter \code{h} is very small (less than 0.001, say), then the approximated observed Fisher information matrix may not be invertible. }

\value{ A list consists of the maximum likelihood estimator, approximated asymptotic standard error, upper, and lower bounds of \eqn{100(1-\alpha)\%} asymptotic confidence interval for parameters of the finite mixture model. }

\references{
K. E. Basford, D. R. Greenway, G. J. McLachlan, and D. Peel, (1997). Standard errors of fitted means under normal mixture, \emph{Computational Statistics}, 12, 1-17.

S. X. Lee and G. J. McLachlan, (2018). EMMIXcskew: An R package for the fitting of a mixture of canonical fundamental skew t-distributions, \emph{Journal of Statistical Software}, 83(3), 1-32, \doi{10.18637/jss.v083.i03}.

M. Teimouri, (2022). Finite mixture of skewed sub-Gaussian stable distributions, \url{https://arxiv.org/abs/2205.14067}.

C. Tortora, R. P. Browne, A. ElSherbiny, B. C. Franczak, and P. D. McNicholas, (2021). Model-based clustering, classification, and discriminant analysis using the generalized hyperbolic distribution: MixGHD R package. \emph{Journal of Statistical Software}, 98(3), 1-24, \doi{10.18637/jss.v098.i03}.
}
\author{Mahdi Teimouri}
\examples{
\donttest{
# Example 1: Approximating the asymptotic standard error and 95 percent confidence interval
#            for the parameters of fitted three-component normal mixture model to iris data.
      Y <- as.matrix( iris[, 1:4] )
colnames(Y) <- NULL
rownames(Y) <- NULL
      G <- 3
 weight <- c( 0.334, 0.300, 0.366         )
    mu1 <- c( 5.0060, 3.428, 1.462, 0.246 )
    mu2 <- c( 5.9150, 2.777, 4.204, 1.298 )
    mu3 <- c( 6.5468, 2.949, 5.482, 1.985 )
 sigma1 <- matrix( c( 0.133, 0.109, 0.019, 0.011, 0.109, 0.154, 0.012, 0.010,
                      0.019, 0.012, 0.028, 0.005, 0.011, 0.010, 0.005, 0.010 ), nrow = 4 , ncol = 4)
 sigma2 <- matrix( c( 0.225, 0.076, 0.146, 0.043, 0.076, 0.080, 0.073, 0.034,
                      0.146, 0.073, 0.166, 0.049, 0.043, 0.034, 0.049, 0.033 ), nrow = 4 , ncol = 4)
 sigma3 <- matrix( c( 0.429, 0.107, 0.334, 0.065, 0.107, 0.115, 0.089, 0.061,
                      0.334, 0.089, 0.364, 0.087, 0.065, 0.061, 0.087, 0.086 ), nrow = 4 , ncol = 4)
     mu <- list(    mu1,    mu2,    mu3 )
  sigma <- list( sigma1, sigma2, sigma3 )
  sigma <- list( sigma1, sigma2, sigma3 )
 lambda <- list( rep(0, 4), rep(0, 4), rep(0, 4) )
   out1 <- sefm( Y, G, weight, model = "restricted", mu, sigma, lambda, family = "constant",
    skewness = "FALSE")
# Example 2: Approximating the asymptotic standard error and 95 percent confidence interval
#            for the parameters of fitted two-component restricted skew t mixture model to
#            AIS data.
      data( AIS )
      Y <- as.matrix( AIS[, 2:3] )
      G <- 2
 weight <- c(  0.5075,  0.4925 )
    mu1 <- c( 19.9827, 17.8882 )
    mu2 <- c( 21.7268,  5.7518 )
 sigma1 <- matrix( c(3.4915, 8.3941, 8.3941, 28.8113 ), nrow = 2, ncol = 2 )
 sigma2 <- matrix( c(2.2979, 0.0622, 0.0622,  0.0120 ), nrow = 2, ncol = 2 )
lambda1 <- ( c( 2.5186, -0.2898 ) )
lambda2 <- ( c( 2.1681,  3.5518 ) )
 theta1 <- c( 68.3088 )
 theta2 <- c(  3.8159 )
     mu <- list(     mu1,     mu2 )
  sigma <- list(  sigma1,  sigma2 )
 lambda <- list( lambda1, lambda2 )
  theta <- list(  theta1,  theta2 )
  param <- c( "nu" )
    PDF <- quote( (nu/2)^(nu/2)*w^(-nu/2 - 1)/gamma(nu/2)*exp( -nu/(w*2) ) )
  tick  <- c( 1, 1 )
   out2 <- sefm( Y, G, weight, model = "restricted", mu, sigma, lambda, family = "igamma",
            skewness = "TRUE", param, theta, tick, h = 0.001, N = 3000, level = 0.05, PDF )
# Example 3: Approximating the asymptotic standard error and 95 percent confidence interval
#            for the parameters of fitted two-component restricted skew sub-Gaussian
#            alpha-stable mixture model to bankruptcy data.
      data( bankruptcy )
      Y <- as.matrix( bankruptcy[, 2:3] );  colnames(Y) <- NULL; rownames(Y) <- NULL
      G <- 2
 weight <- c(  0.553,  0.447 )
    mu1 <- c( -3.649, -0.085 )
    mu2 <- c( 40.635, 19.042 )
 sigma1 <- matrix( c(1427.071, -155.356, -155.356, 180.991 ), nrow = 2, ncol = 2 )
 sigma2 <- matrix( c( 213.938,    9.256,    9.256,  74.639 ), nrow = 2, ncol = 2 )
lambda1 <- c( -41.437, -21.750 )
lambda2 <- c(  -3.666,  -1.964 )
 theta1 <- c( 1.506 )
 theta2 <- c( 1.879 )
     mu <- list(     mu1,     mu2 )
  sigma <- list(  sigma1,  sigma2 )
 lambda <- list( lambda1, lambda2 )
  theta <- list(  theta1,  theta2 )
  param <- c( "alpha" )
  tick  <- c( 1 )
   out3 <- sefm( Y, G, weight, model = "restricted", mu, sigma, lambda, family = "pstable",
            skewness = "TRUE", param, theta, tick, h = 0.01, N = 3000, level = 0.05 )
# Example 4: Approximating the asymptotic standard error and 95 percent confidence interval
#            for the parameters of fitted two-component restricted generalized inverse-Gaussian
#            mixture model to AIS data.
      data( wheat )
      Y <- as.matrix( wheat[, 1:7] ); colnames(Y) <- NULL; rownames(Y) <- NULL
      G <- 3
 weight <- c( 0.325, 0.341, 0.334 )
    mu1 <- c( 18.8329, 16.2235, 0.9001, 6.0826, 3.8170, 1.6604, 6.0260 )
    mu2 <- c( 11.5607, 13.1160, 0.8446, 5.1873, 2.7685, 4.9884, 5.2203 )
    mu3 <- c( 13.8071, 14.0720, 0.8782, 5.5016, 3.1513, 0.6575, 4.9111 )
lambda1 <- diag( c( 0.1308, 0.2566,-0.0243, 0.2625,-0.1259, 3.3111, 0.1057) )
lambda2 <- diag( c( 0.7745, 0.3084, 0.0142, 0.0774, 0.1989,-1.0591,-0.2792) )
lambda3 <- diag( c( 2.0956, 0.9718, 0.0042, 0.2137, 0.2957, 3.9484, 0.6209) )
 theta1 <- c( -3.3387, 4.2822 )
 theta2 <- c( -3.6299, 4.5249 )
 theta3 <- c( -3.9131, 5.8562 )
 sigma1 <- matrix( c(
 1.2936219, 0.5841467,-0.0027135, 0.2395983, 0.1271193, 0.2263583, 0.2105204,
 0.5841467, 0.2952009,-0.0045937, 0.1345133, 0.0392849, 0.0486487, 0.1222547,
-0.0027135,-0.0045937, 0.0003672,-0.0033093, 0.0016788, 0.0056345,-0.0033742,
 0.2395983, 0.1345133,-0.0033093, 0.0781141, 0.0069283,-0.0500718, 0.0747912,
 0.1271193, 0.0392849, 0.0016788, 0.0069283, 0.0266365, 0.0955757, 0.0002497,
 0.2263583, 0.0486487, 0.0056345,-0.0500718, 0.0955757, 1.9202036,-0.0455763,
 0.2105204, 0.1222547,-0.0033742, 0.0747912, 0.0002497,-0.0455763, 0.0893237 ), nrow = 7, ncol = 7 )
 sigma2 <- matrix( c(
 0.9969975, 0.4403820, 0.0144607, 0.1139573, 0.1639597,-0.2216050, 0.0499885,
 0.4403820, 0.2360065, 0.0010769, 0.0817149, 0.0525057,-0.0320012, 0.0606147,
 0.0144607, 0.0010769, 0.0008914,-0.0023864, 0.0049263,-0.0122188,-0.0042375,
 0.1139573, 0.0817149,-0.0023864, 0.0416206, 0.0030268, 0.0490919, 0.0407972,
 0.1639597, 0.0525057, 0.0049263, 0.0030268, 0.0379771,-0.0384626,-0.0095661,
-0.2216050,-0.0320012,-0.0122188, 0.0490919,-0.0384626, 4.0868766, 0.1459766,
 0.0499885, 0.0606147,-0.0042375, 0.0407972,-0.0095661, 0.1459766, 0.0661900 ), nrow = 7, ncol = 7 )
 sigma3 <- matrix( c(
 1.1245716, 0.5527725,-0.0005064, 0.2083688, 0.1190222,-0.4491047, 0.2494994,
 0.5527725, 0.3001219,-0.0036794, 0.1295874, 0.0419470,-0.1926131, 0.1586538,
-0.0005064,-0.0036794, 0.0004159,-0.0034247, 0.0019652,-0.0026687,-0.0044963,
 0.2083688, 0.1295874,-0.0034247, 0.0715283, 0.0055925,-0.0238820, 0.0867129,
 0.1190222, 0.0419470, 0.0019652, 0.0055925, 0.0243991,-0.0715797, 0.0026836,
-0.4491047,-0.1926131,-0.0026687,-0.0238820,-0.0715797, 1.5501246,-0.0048728,
 0.2494994, 0.1586538,-0.0044963, 0.0867129, 0.0026836,-0.0048728, 0.1509183 ), nrow = 7, ncol = 7 )
    mu <- list( mu1, mu2, mu3 )
 sigma <- list( sigma1 , sigma2, sigma3 )
lambda <- list( lambda1, lambda2, lambda3 )
 theta <- list( theta1 , theta2, theta3 )
  tick <- c( 1, 1, 0 )
 param <- c( "a", "b" )
   PDF <- quote( 1/( 2*besselK( b, a ) )*w^(a - 1)*exp( -b/2*(1/w + w) ) )
  out4 <- sefm( Y, G, weight, model = "unrestricted", mu, sigma, lambda, family = "gigaussian",
            skewness = "TRUE", param, theta, tick, h = 0.001, N = 3000, level = 0.05, PDF )
}
}
