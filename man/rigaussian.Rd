\name{rigaussian}
\alias{rigaussian}
\title{Simulating from inversse Gaussian random variable.}
\description{Using method of Michael and Schucany (1976), we can generate from inversse Gaussian random variable. The density function of an inversse Gaussian distribution is given by
\deqn{ f_W(w|{\bold{\theta}}) =\sqrt{\frac{\beta}{2 \pi w^3}}\exp\biggl\{-\frac{\beta(w - \mu)^2}{2\mu^2 w}\biggr\},} where \eqn{w>0} and \eqn{{\bold{\theta}}=(\mu, \beta)^{\top}}. Herein \eqn{\mu>0} is the mean and \eqn{\beta> 0} are the first (mean) and second (shape) parameter of this family, respectively.}
\usage{rigaussian(n, alpha, beta)}
\arguments{
\item{n}{size of required samples.}
\item{mu}{tail mean parameter.}
\item{beta}{ shape parameter.}
}
%\details{}
\value{ simulated realizations of size \eqn{n} from inversse Gaussian random variable. }
\references{
J. R. Michael and Schucany, (1976). Generating Random Variates Using Transformations with Multiple Roots, \emph{The American Statistician}, 30(2), 88-90, \doi{10.1080/00031305.1976.10479147}.
}
\author{Mahdi Teimouri}
\examples{
\donttest{
    n <- 100
alpha <- 4
 beta <- 2
rigaussian(n, alpha, beta)
}
}
\keyword{internal}
