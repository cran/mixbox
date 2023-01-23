rptstable <- function(n, alpha, beta)
{
  if( alpha <= 0 | alpha >2 )  { stop (gettextf("parameter alpha must be in (0, 2].")) }
  if( beta <= 0             )  { stop (gettextf("parameter beta  must be positive.")) }
  a  <- alpha/2
  b0 <- beta
  b1 <- b0/a
  b2 <- a^(-a)*(1 - a)^(-(1 - a))
  b3 <- 1/sqrt(b1*a*(1 - a))
  b4 <- rep(NA, n)
  y  <- rep(NA, n)
  i  <- 1
  while(i <= n)
  {
    if(b3 >= sqrt(2 * pi))
    {
      u0 <- runif(1, 0, pi)
      u1 <- runif(1)
      u2 <- ( sin(u0)/( ( sin(a*u0) )^(a)*( sin((1 - a)*u0))^(1 - a)))^b1
      if ( (u1*b2^b1) < u2 )
      {
        y[i] <- u0
        i <- i + 1
      }
    }else{
      N  <- rnorm(1)
      u1 <- runif(1)
      u2 <- ( sin(b3*abs(N) )/(( sin(a*b3*abs(N)) )^(a)*( sin((1 - a)*b3*abs(N)) )^(1 - a)) )^b1
      if ( ( (b3*abs(N) ) < pi ) && (u1*b2^b1*exp( -N^2/2 )) < u2 )
      {
        y[i] <- b3*abs(N)
        i <- i + 1
      }
    }
  }
  b4  <- ( ( (sin(a*y))^(a)*(sin((1 - a)*y))^(1 - a) )/sin(y) )^( 1/(1 - a) )
  out <- ( b4/rgamma(n, (1 + 2*b0*(1 - a)/(2*a)), 1))^(2*(1 - a)/(2*a) )
  return( out )
}
