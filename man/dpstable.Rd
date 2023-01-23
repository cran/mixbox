\name{dpstable}
\alias{dpstable}
\title{Monte Carlo approximation for density function of positive alpha-stable distribution.}
\description{The density function \eqn{f_{P}(p|\alpha)}, of positive \eqn{\alpha}-stable distribution is given by (Kanter, 1975):
\deqn{
 f_{P}(p|\alpha)=\frac{1}{\pi}\int_{0}^{\pi}{\frac{\alpha}{2-\alpha}}a(\theta) p^{-\frac{\alpha}{2-\alpha}-1}a(\theta) \exp\Bigl\{-p^{-\frac{\alpha}{2-\alpha}}a(\theta)\Bigr\}d\theta, }
where \eqn{0<\alpha \leq 2} is tail thickness parameter or index of stability and
\deqn{
a(\theta)=\frac{\sin\Bigl(\bigl(1-\frac{\alpha}{2}\bigr)\theta\Bigr)\Bigl[\sin \bigl(\frac{\alpha \theta}{2}\bigr)\Bigr]^{\frac{\alpha}{2-\alpha}}}{[\sin(\theta)]^{\frac{2}{2-\alpha}},}
}
for \eqn{0<\theta < \pi}. We use the Monte Carlo method for approximating \eqn{f_{P}(p|\alpha)}.
}
\usage{dpstable(x, param)}
\arguments{
\item{x}{point at which density value is desired. }
\item{param}{tail index parameter. }
}

\value{
The density function of positive \eqn{\alpha}-stable distribution at point \eqn{x}.
}
\references{
M. Kanter, (1975). Stable densities under change of scale and total variation inequalities, \emph{Annals of Probability}, 3(4), 697-707.
}
\author{Mahdi Teimouri}
\examples{
\donttest{
x <- 2
param <- 1.5
dpstable(x, param)
}
}
\keyword{internal}
