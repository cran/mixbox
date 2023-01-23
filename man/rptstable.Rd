\name{rptstable}
\alias{rptstable}
\title{Simulating from polynomially tilted \eqn{\alpha}-stable (ptstable) random variable.}
\description{Using method of Devroye (2009), we can generate from ptstable random variable. The density function of a ptstable distribution is given by
\deqn{
f_{T}(t|{\bold{\theta}})=\frac{\Gamma(1+\beta)}{\Gamma\Bigl(1+\frac{\beta}{\alpha}\Bigr)}t^{-\beta}f_{P}(t|\alpha),
}
where \eqn{{\bold{\theta}}=(\alpha,\beta)^{\top}} in which \eqn{0<\alpha \leq 2} is tail thickness parameter or index of stability and \eqn{\beta> 0} is tilting parameter. We note that \eqn{f_{P}(t|\alpha)} is the density function of a positive \eqn{\alpha}-stable distribution that has an integral representation (Kanter, 1975):
\deqn{ f_{P}(t|\alpha)=\frac{1}{\pi}\int_{0}^{\pi}{\frac{\alpha}{2-\alpha}}a(\theta) t^{-\frac{\alpha}{2-\alpha}-1}a(\theta) \exp\Bigl\{-t^{-\frac{\alpha}{2-\alpha}}a(\theta)\Bigr\}d\theta, }
for
\deqn{
a(\theta)=\frac{\sin\Bigl(\bigl(1-\frac{\alpha}{2}\bigr)\theta\Bigr)\Bigl[\sin \bigl(\frac{\alpha \theta}{2}\bigr)\Bigr]^{\frac{\alpha}{2-\alpha}}}{[\sin(\theta)]^{\frac{2}{2-\alpha}}}.
}
}
\usage{rptstable(n, alpha, beta)}
\arguments{
\item{n}{size of required samples.}
\item{alpha}{tail thickness parameter.}
\item{beta}{ tilting parameter.}
}
%\details{}
\value{ simulated realizations of size \eqn{n} from ptstable distribution. }
\references{
M. Kanter, (1975). Stable densities under change of scale and total variation inequalities, \emph{Annals of Probability}, 3(4), 697-707.

L. Devroye, (2009). Random variate generation for exponentially and polynomially tilted stable distributions, \emph{ACM Transactions on Modeling and Computer Simulation}, 19(4), \doi{10.1145/1596519.1596523}.
}
\author{Mahdi Teimouri}
\examples{
\donttest{
    n <- 100
alpha <- 1.4
beta  <- 0.5
rptstable(n, alpha, beta)
}
}
\keyword{internal}
