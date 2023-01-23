dpstable <- function(x, param)
{
  n   <- length( x )
  out <- rep( NA, n )
 # n.k <- 90
 # m   <- 1500
 # k   <- seq( 1, n.k )
 # L   <- 5 + 2*( exp( lgamma( param*n.k/2 + param/2 + 1 ) + lgamma( param*n.k/2 + param/2 + 3/2 ) -                      lgamma( param*n.k/2 + 1 ) - lgamma( param*n.k/2 + 3/2 ) )/( (n.k + 1) ) )^( 2/param )
 # for(j in 1: n)  {    if(x[j] > L)    {      out[j] <- 1/(pi)*sum( (-1)^(k - 1)*exp(lgamma( param*k/2 + 1 ) - lgamma(k + 1) )*sin(k*pi*param/2)*x[j]^( -param*k/2 - 1 )  )    }else{
     # theta  <- runif(m, 0, pi);      out[j] <- mean( sapply( 1:m, function(i)                                      {                                        param/( 2 - param )*x[j]^( -param/( 2 - param ) - 1)*sin( (1 - param/2)*theta[i] )*                                        ( sin(param*theta[i]/2) )^( param/(2 - param ) )/( sin(theta[i]) )^( 2/(2 - param ) )*                                        exp( -x[j]^( -param/( 2 - param ) )*sin( (1 - param/2 )*theta[i] )*( sin(param*theta[i]/2) )^                                        ( param/(2 - param) )/( sin(theta[i]) )^( 2/(2 - param ) )  ) }       )                    )
     # out[j] <- dstable(x[j],param/2,1,cos(pi*param/4)^(2/param),0,1)

 #   }
 # }
  out <- suppressWarnings( dstable(x, param/2, 1, cos(pi*param/4)^(2/param), 0, 1) )
  return( out )
}
