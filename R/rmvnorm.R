rmvnorm <- function(n, Mu, Sigma)
{
  Dim <- length( Mu )
    Y <- matrix(NA, nrow = n, ncol = Dim)
    L <- t( chol( Sigma ) )
      for(i in 1:n)
      {
        Y[i, ] <- L %*% rnorm( Dim ) + Mu
      }
    return( Y )
}
