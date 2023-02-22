arrange_sigma <- function(x)
{
  Dim <- length( x[1, ] )
  y <- rep( 0, Dim*(Dim + 1)/2 )
  index <- matrix(0, nrow = Dim*(Dim + 1)/2, ncol = 2)
  k <- 1
  for(i in 1:Dim)
  {
    for(j in i:Dim)
    {
      y[k] <- x[i, j]
      index[k, ] <- c(i, j)
      k <- k + 1
    }
  }
return( list(y = y, index = index) )
}
