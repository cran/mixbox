rmix <- function(n, G, weight, model = "restricted", mu, sigma, lambda, family = "constant", theta = NULL)
{
  if( G < 1 | G != round(G) ) stop( "Number of components must be integer." )
    if( G == 1 & length( weight ) > 1 ) stop( "Length of the mixing proportions must be one." )
      if( sum( weight ) != 1 || any( weight < 0 ) ) stop( "Elements of mixing proportion must be positive and must sum to 1." )
        if( length( weight ) != G ) stop( "Length of mixing proportions and number of components G must be equal." )
          if( length( mu     ) != G ) stop( "Length of the location parameter and number of components G must be equal." )
            if( length( sigma  ) != G ) stop( "Length of the dispersion matrix and number of components G must be equal." )
              if( length( lambda ) != G ) stop( "Length of the skewness parameter and number of components G must be equal." )
                if( length( theta  ) != G ) stop("The length of the ML estimators of mixing distrinution and number of components G must be equal.")
                  if(family != "constant" & family != "bs" & family != "burriii" & family != "chisq" & family != "exp" & family != "f" &
                    family != "gamma" & family != "gigaussian"  & family != "igamma" & family != "lidley" & family != "loglog" &
                    family != "lognorm" & family != "lomax" & family != "pstable" & family != "ptstable" & family != "rayleigh" & family != "weibull" )
                    stop( "Mixing distribution not implemented or misspelled. Please check the manual for guidelines." )
                  if( all( model != c( "canonical", "restricted", "unrestricted" ) ) ) stop( "model's name must be either canonical, restricted, or unrestricted." )
                if( family == "constant" ) theta <- rep( list(200), G )
              if( family != "constant" & is.null( theta ) ) stop( "ML estimator of the mixing distribution parameters must be given." )
            for(g in 1:G)
            {
              if( model == "restricted" & is.numeric( dim( lambda[[g]] ) )  ) stop("The skewness parameter must be a vector of appropriate size.")
            }
          for(g in 1:G)
          {
            if( model == "unrestricted" & !is.numeric( dim( lambda[[g]] ) )  ) stop("The skewness parameter must be a diagonal matrix of appropriate dimension.")
          }
        if( model == "unrestricted" &   any( dim( lambda[[1]]  ) != dim( sigma[[1]] ) ) ) stop("skewness and dispersion matrices muust be of the same size.")
      if( model == "unrestricted" |   model == "canonical"   ) Q <- length( lambda[[g]][1, ] )
    for( g in 1:G )
    {
    if( any( eigen( sigma[[g]] )$values <= 0 ) ) stop("dispersion matrix must be positive definite.")
    }
  Dim  <- length( mu[[1]] )
  W    <- rep( NA, n )
  Z0   <- rep( NA, n )
  n_g  <- rep( NA, G )
  cn_g <- rep( NA, G )
    Y  <- matrix( NA, nrow = n, ncol = Dim )
  Z1   <- matrix( NA, nrow = n, ncol = Dim )
  sim <- rmultinom(n, 1, weight)
  n_g  <- apply( sim, 1, sum )
  label <- rep( seq(1, G), n_g )
  cn_g <- cumsum( n_g )
  for(g in 1:G)
  {
    if( family == "constant")
    {
      Dim_theta <- 1
      param     <-  c( "nu" )
      W         <- rep( 1, n_g[g] )
    }
    if( family != "constant" )
    {
      if(family == "bs")
      {
        z0 <- theta[1]*rnorm(n_g[g])/2
        W  <- theta[2]*( 1 + 2*z0^2 + 2*z0*sqrt(1 + z0^2) )
      }
      if(family == "burriii")
      {
        W <- exp( -log( exp(- log( runif( n_g[g] ) )/theta[[g]][1] ) - 1  )/theta[[g]][2] )
      }
      if(family == "chisq")
      {
        W <- rchisq( n_g[g], df = theta[[g]][1] )
      }
      if(family == "f")
      {
        W <- rf( n_g[g], df1 = theta[[g]][1], df2 = theta[[g]][2] )
      }
      if(family == "gamma")
      {
        W <- rgamma( n_g[g], shape = theta[[g]][1], rate = theta[[g]][2] )
      }
      if(family == "igaussian")
      {
        W <- rigaussian( n_g[g], alpha = theta[[g]][1], beta = theta[[g]][2] )
      }
      if(family == "igamma")
      {
        W <- 1/rgamma( n_g[g], shape = theta[[g]][1], rate = theta[[g]][2] )
      }
      if(family == "gigaussian")
      {
        W <- rgig( n_g[g], lambda = theta[[g]][1], chi = theta[[g]][2], psi = theta[[g]][3] )
      }
      if(family == "lidley")
      {
        x <- rep(NA, n_g[g])
        x <- rbinom( n_g[g], 1, theta[[g]][1]/(1 + theta[[g]][1]) )
        W <- x*rgamma( n_g[g], shape = 1, rate = theta[[g]][1] ) + (1 - x)*rgamma( n_g[g], shape = 2, rate = theta[[g]][1] )
      }
      if(family == "loglog")
      {
        u <- runif( n_g[g] )
        W <- theta[[g]][2]*( u/(1 - u) )^(1/theta[[g]][1])
      }
      if(family == "lognorm")
      {
        W <- rlnorm( n_g[g], meanlog = theta[[g]][1], sdlog = theta[[g]][2] )
      }
      if(family == "lomax")
      {
        W <- 1/theta[[g]][2]*( (1 - runif( n_g[g] ) )^(-1/theta[[g]][1]) - 1 )
      }
      if(family == "pstable")
      {
        w0     <- -log( runif(n_g[g]) )
        theta0 <- pi*(runif(n_g[g]) - 1/2)
        sigma0 <- (cos(pi*theta[[g]][1]/4))^(2/theta[[g]][1])
        r0     <- sin(theta[[g]][1]/2*(pi/2 + theta0))/( cos(theta[[g]][1]*pi/4)*cos(theta0) )^(2/theta[[g]][1])*
          ( cos(theta[[g]][1]*pi/4 + (theta[[g]][1]/2 - 1)*theta0)/w0 )^( (2 - theta[[g]][1])/theta[[g]][1] )
        W <- ( r0 - tan(pi*theta[[g]][1]/4) )*sigma0 + sigma0*tan(pi*theta[[g]][1]/4)
      }
      if(family == "ptstable")
      {
        W <- 1/( 2^(1 - 2/theta[[g]][1])*rptstable(n_g[g], theta[[g]][1], Dim/2) )
      }
      if(family == "rayleigh")
      {
        W <- rweibull( n_g[g], shape = 2, scale = theta[[g]][1] )
      }
      if(family == "weibull")
      {
        W <- rweibull( n_g[g], shape = theta[[g]][1], scale = theta[[g]][2] )
      }
    }
    Z1  <- rmvnorm( n = n_g[g], Mu = rep(0, Dim), Sigma = sigma[[g]] )
    Mu     <- as.vector( mu[[g]] )
    Sigma  <- as.matrix( sigma[[g]] )
    Lambda <- as.matrix( lambda[[g]] )
    if ( any( model == c("canonical", "unrestricted") ) )
    {
      Z0  <- matrix( qnorm( ( runif(n_g[g]*Q) + 1 )/2 ), nrow = n_g[g], ncol = Q, byrow = TRUE )
      Y[ ( 1 + sign( g - 1 )*cn_g[ g - sign(g - 1) ] ):cn_g[g], ] <- t( Mu + sapply(1:n_g[g], function(i) sqrt( W[i] )*( Lambda%*%Z0[i, ] + Z1[i, ] ) ) )
    }else{
      Z0  <- qnorm( ( runif(n_g[g]) + 1 )/2 )
      Y[ ( 1 + sign(g - 1)*cn_g[ g - sign(g - 1) ] ):cn_g[g], ] <- t( Mu + sapply(1:n_g[g], function(i) sqrt(W[i])*( Lambda*Z0[i] + Z1[i, ] ) ) )
    }
  }
  out <- cbind( Y, label )
  colnames( out ) <- NULL
  return( out )
}
