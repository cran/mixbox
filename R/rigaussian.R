rigaussian <- function(n, alpha, beta)
{
      y <- rep(NA, n)
     z2 <-rnorm(n)^2
  ratio <- alpha + alpha^2*z2/(2*beta) - alpha/(2*beta)*sqrt( 4*alpha*beta*z2 + (alpha*z2)^2 )
      y <- alpha^2/ratio
      u <- runif(n)
      j <- which( ratio < alpha*u/(1 - u) )
   y[j] <- ratio[j]
return(y)
}
