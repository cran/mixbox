arrange_lambda<- function(x)
{
  Dim <- dim( x )[1]
  Q   <- dim( x )[2]
  y   <- rep( 0, Dim*Q )
  index <- matrix(0, nrow = Dim*Q, ncol = 2)
  k <- 1
  for(i in 1:Dim)
  {
    for(j in 1:Q)
    {
      y[k] <- x[i, j]
      index[k, ] <- c(i, j)
      k <- k + 1
    }
  }
return( list(y = y, index = index) )
}
