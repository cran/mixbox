\name{dptstable}
\alias{dptstable}
\title{Monte Carlo approximation for density function of polynomially tilted alpha-stable distribution.}
\description{The density function \eqn{f_{T}(t|\alpha, \beta)}, of polynomially tilted \eqn{\alpha}-stable distribution is given by (Devroye, 2009):
\deqn{f_{T}(t | \alpha, \beta)=\frac{\Gamma(1+\beta)}{\Gamma\Bigl(1+\frac{\beta}{\alpha}\Bigr)}t^{-\beta}f_{P}(t|\alpha),}
where \eqn{0<\alpha \leq 2} is tail thickness parameter or index of stability and \eqn{\beta> 0} is tilting parameter. We note that \eqn{f_{P}(t|\alpha)} is the density function of a positive \eqn{\alpha}-stable distribution that has an integral representation (Kanter, 1975):
\deqn{ f_{P}(t|\alpha)=\frac{1}{\pi}\int_{0}^{\pi}{\frac{\alpha}{2-\alpha}}a(\theta) t^{-\frac{\alpha}{2-\alpha}-1}a(\theta) \exp\Bigl\{-t^{-\frac{\alpha}{2-\alpha}}a(\theta)\Bigr\}d\theta, }
where
\deqn{
a(\theta)=\frac{\sin\Bigl(\bigl(1-\frac{\alpha}{2}\bigr)\theta\Bigr)\Bigl[\sin \bigl(\frac{\alpha \theta}{2}\bigr)\Bigr]^{\frac{\alpha}{2-\alpha}}}{[\sin(\theta)]^{\frac{2}{2-\alpha}}},
}
for \eqn{0 < \theta < \pi}.
}
\usage{dptstable(x, param, Dim)}
\arguments{
\item{x}{point at which density value is desired.}
\item{param}{tail thickness parameter.}
\item{Dim}{tilting parameter.}
}
%\details{Herein, we consider \eqn{\beta = d/2} where \eqn{d} denotes the dimension of random vector \eqn{\bold{Y}=}. }
\value{ The density function of polynomially tilted \eqn{\alpha}-stable distribution at point \eqn{x}.}
\references{
M. Kanter, (1975). Stable densities under change of scale and total variation inequalities, \emph{Annals of Probability}, 3(4), 697-707.

L. Devroye, (2009). Random variate generation for exponentially and polynomially tilted stable distributions, \emph{ACM Transactions on Modeling and Computer Simulation}, 19(4), \doi{10.1145/1596519.1596523}.
}
\author{Mahdi Teimouri}
\examples{
\donttest{
    x <- 2
param <- 1.5
  Dim <- 2
dptstable(x, param, Dim)
}
}
\keyword{internal}
