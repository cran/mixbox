dptstable <- function(x, param, Dim)
{
  n   <- length( x )
  out <- rep( NA, n )
 # n.k <- 80
 # m   <- 500
 # k   <- seq( 1, n.k )
 # L   <- 3 + 2*( exp( lgamma( param[1]*n.k/2 + param[1]/2 + 1 ) + lgamma( param[1]*n.k/2 + param[1]/2 + 3/2 ) -
 #                     lgamma( param[1]*n.k/2 + 1 ) - lgamma( param[1]*n.k/2 + 3/2 ) )/( (n.k + 1) ) )^( 2/param[1] )
 #  for(j in 1: n)
 # {
 #   if(x[j] > L)
 #   {
 #     out[j] <- 1/(pi)*sum( (-1)^(k - 1)*exp(lgamma( param[1]*k/2 + 1 ) - lgamma(k + 1) )*sin(k*pi*param[1]/2)*( 2^( 2/param[1] - 1 )*x[j] )^( -param[1]*k/2 - 1 )  )
 #   }else{
 #     theta  <- runif(m, 0, pi)
 #      out[j] <- mean( sapply( 1:m, function(i)
 #                                     {
 #                                       param[1]/( 2 - param[1] )*x[j]^( -param[1]/( 2 - param[1] ) - 1)*sin( (1 - param[1]/2)*theta[i] )*
 #                                       ( sin(param[1]*theta[i]/2) )^( param[1]/(2 - param[1]) )/( sin(theta[i]) )^( 2/(2 - param[1]) )*
 #                                       exp( -x[j]^( -param[1]/( 2 - param[1] ) )*sin( (1 - param[1]/2)*theta[i] )*( sin(param[1]*theta[i]/2) )^
 #                                       ( param[1]/(2 - param[1]) )/( sin(theta[i]) )^( 2/(2 - param[1]) ) )
 #                                     }
 #                           )
 #                   )
 #   }
 # }
  out <- dstable(x,param/2,1,cos(pi*param/4)^(2/param),0,1)
  2^( Dim/2 - Dim/param[1] )*gamma( 1 + param[1]/2 )/gamma( 1 + Dim/param[1] )*x^(-Dim/2)*out
}
