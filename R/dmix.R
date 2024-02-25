dmix <- function(Y, G, weight, model = "restricted", mu, sigma, lambda, family = "constant", skewness = "FALSE", param = NULL, theta = NULL, tick = NULL, N = 3000, log = "FALSE")
{
if( G < 1 | G != round(G) ) stop( "Number of components must be integer." )
  if( G == 1 & length( weight ) > 1 ) stop( "Length of the mixing proportions must be one." )
    if( sum( weight ) != 1 || any( weight < 0 ) ) stop( "Elements of mixing proportion must be positive and must sum to 1." )
      if( length( weight ) != G ) stop( "Length of mixing proportions and number of components G must be equal." )
        if( length( mu     ) != G ) stop( "Length of the location parameter and number of components G must be equal." )
          if( length( sigma  ) != G ) stop( "Length of the dispersion matrix and number of components G must be equal." )
            if( length( lambda ) != G ) stop( "Length of the skewness parameter and number of components G must be equal." )
              if(family != "constant" & family != "bs" & family != "burriii" & family != "chisq" & family != "exp" &
                family != "f" & family != "gamma" & family != "gigaussian"  & family != "igamma" & family != "lidley" &
                  family != "loglog" & family != "lognorm" & family != "lomax" & family != "pstable" & family != "ptstable" &
                    family != "rayleigh" & family != "weibull" )
                      { stop( "Mixing distribution not implemented or misspelled. Please check the manual for guidelines." ) }
                        if( skewness != TRUE & skewness != FALSE ) stop( "Skewness must be a logical statement either TRUE or FALSE." )
                          if( family != "constant" & is.null( theta ) ) stop( "ML estimator of the mixing distribution parameters must be given." )
                        if( family != "constant" & is.null( tick ) ) stop( "vector tick must be given." )
                      if( family != "constant" & length( theta ) != G ) stop( "Length of the ML estimators of mixing distribution and number of
	                  components G must be equal." )
                  if( length( param ) != length( theta[[1]] ) ) stop( "Length of the parameter vector of mixing distribution and associated
	              MLEs must be equal." )
              if(family != "constant")
              {
              if( any( tick < 0 ) || any( tick > 1 ) || sum( tick ) == 0 || ( sum( tick )%%1 != 0 & exp( prod( tick ) ) != exp(1) ) ) stop( "Elements
                                                                                                                of vector tick are either 0 or 1." )
              }
          if( family != "constant" & length(tick) < length(param) ) stop( "Length of vector tick cannot exceed the length of param." )
        if( all( model != c("canonical", "restricted", "unrestricted") ) ) stop( "model must be canonical, restricted, or unrestricted." )
      if(model == "restricted") Q <- length( lambda[[1]] )
    if(model == "canonical" | model == "unrestricted" ) Q <- dim( lambda[[1]] )[2]
  Dim <- length( mu[[1]] )
  n   <- ifelse( is.null( dim(Y) ), 1, length( Y[, 1] ) )
  if( n == 1 ) Y <- t( as.matrix( c(Y), nrow = 1, ncol = Dim ) )
  if(model == "restricted")
  {
    s1 <- rep(1, G)
    if( skewness == "TRUE" )
    {
      for(g in 1:G)
      {
        if( all( lambda[[g]] == 0 ) ) s1[g] <- 0
      }
    }else{
      lambda <- vector("list", G)
      for(g in 1:G)
      {
        lambda[[g]] <- rep( 0, Dim )
      }
    }
    if( sum( s1 ) == 0 ) stop( "Skewness vector must be non-zero." )
  }
  if(model == "canonical")
  {
    s2 <- rep(1, G)
    s3 <- rep(1, G)
    if( skewness == "TRUE" )
    {
      for(g in 1:G)
      {
        for(i in 1:Dim)
        {
          if( all( lambda[[g]][i, ] == 0 ) ) s2[g] <- 0
        }
        for(j in 1:Q)
        {
          if( all( lambda[[g]][ , j] == 0 ) ) s3[g] <- 0
        }
      }
    }else{
      lambda <- vector("list", G)
      for(g in 1:G)
      {
        lambda[[g]] <- matrix( 0, nrow = Dim, ncol = Q )
      }
    }
    if( sum( s2 ) == 0 | sum( s3 ) == 0 ) stop( "Skewness vector must be non-zero." )
  }
  P1 <- C0 <- delta <- rep( 0, G )
  dy <- m  <- f_y   <- matrix( NA, nrow = n, ncol = G )
  W  <- matrix( NA, nrow = N, ncol = G )
  if( model != "restricted" )
  {
    T0 <- array( NA, c( N, Q, G ) )
    T1 <- array( NA, c( N, Q, G ) )
  }
  for (g in 1:G)
  {
    if(family == "burriii")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- exp( -log( exp(- log( runif( N ) )/theta[[g]][1] ) - 1 )/theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- exp( -log( exp(- log( runif( N ) )/theta[[g]][1] ) - 1 )/theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- exp( -log( exp(- log( runif( N ) )/1             ) - 1 )/theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <- exp( -log( exp(- log( runif( N ) )/theta[[g]][1] ) - 1 )/1             )
    }
    if(family == "bs")
    {
      if( length(param) == 2 & sum( tick ) == 2 )
      {
        z0 <- theta[[g]][1]*rnorm(N)/2
        W[, g] <- theta[[g]][2]*( 1 + 2*z0^2 + 2*z0*sqrt(1 + z0^2) )
      }
      if( length(param) == 1 & sum( tick ) == 2 )
      {
        z0 <- theta[[g]][1]*rnorm(N)/2
        W[, g] <- theta[[1]][2]*( 1 + 2*z0^2 + 2*z0*sqrt(1 + z0^2) )
      }
      if( length(param) == 1 & tick[1] == 0 )
      {
        z0 <- 1*rnorm(N)/2
        W[, g] <- theta[[1]][1]*( 1 + 2*z0^2 + 2*z0*sqrt(1 + z0^2) )
      }
      if( length(param) == 1 & tick[2] == 0 )
      {
        z0 <- theta[[g]][1]*rnorm(N)/2
        W[, g] <- 1*( 1 + 2*z0^2 + 2*z0*sqrt(1 + z0^2) )
      }
    }
    if(family == "chisq") W[, g] <- rchisq(N, df = theta[[g]][1])
    if(family == "f")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- rf( N, df1 = theta[[g]][1], df2 = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- rf( N, df1 = theta[[g]][1], df2 = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- rf( N, df1 = 1            , df2 = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <- rf( N, df1 = theta[[g]][1], df2 = 1             )
    }
    if(family == "gamma")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- rgamma( N, shape = theta[[g]][1], rate = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- rgamma( N, shape = theta[[g]][1], rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- rgamma( N, shape = 1            , rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <- rgamma( N, shape = theta[[g]][1], rate = 1             )
    }
    if(family == "igamma")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- 1/rgamma( N, shape = theta[[g]][1], rate = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- 1/rgamma( N, shape = theta[[g]][1], rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- 1/rgamma( N, shape = 1            , rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <- 1/rgamma( N, shape = theta[[g]][1], rate = 1             )
    }
    if(family == "gigaussian")
    {
      if( length(param) == 3 & sum( tick ) == 3 ) W[, g] <- rgig( N, lambda = theta[[g]][1], chi = theta[[g]][2], psi = theta[[g]][3] )
      if( length(param) == 2 & sum( tick ) == 3 ) W[, g] <- rgig( N, lambda = theta[[g]][1], chi = theta[[g]][2], psi = theta[[g]][2] )
      if( length(param) == 2 &     tick[1] == 0 ) W[, g] <- rgig( N, lambda = 1, chi = theta[[g]][1], psi = theta[[g]][2] )
      if( length(param) == 2 &     tick[2] == 0 ) W[, g] <- rgig( N, lambda = theta[[g]][1], chi = 1, psi = theta[[g]][2] )
      if( length(param) == 2 &     tick[3] == 0 ) W[, g] <- rgig( N, lambda = theta[[g]][1], chi = theta[[g]][2], psi = 1 )
    }
    if(family == "igaussian")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- rigaussian( N, alpha = theta[[g]][1], beta = theta[[g]][2] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- rigaussian( N, alpha = 1            , beta = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <- rigaussian( N, alpha = theta[[g]][1], beta = 1             )
    }
    if(family == "lidley")
    {
      x <- rep( NA, N )
      x <- rbinom( N, 1, theta[[g]][1]/(1 + theta[[g]][1]) )
      W[, g] <- x*rgamma( N, shape = 1, rate = theta[[g]][1] ) + (1 - x)*rgamma( N, shape = 2, rate = theta[[g]][1] )
    }
    if(family == "loglog")
    {
      u <- runif( N )
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- theta[[g]][2]*( u/(1 - u) )^(1/theta[[g]][1] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- theta[[g]][1]*( u/(1 - u) )^(1/theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- theta[[g]][2]*( u/(1 - u) )^(1/1             )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <-             1*( u/(1 - u) )^(1/theta[[g]][1] )
    }
    if(family == "lognorm")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- rlnorm( N, meanlog = theta[[g]][1], sdlog = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- rlnorm( N, meanlog = theta[[g]][1], sdlog = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- rlnorm( N, meanlog = 0            , sdlog = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <- rlnorm( N, meanlog = theta[[g]][1], sdlog = 1             )
    }
    if(family == "lomax")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- 1/theta[[g]][2]*( (1 - runif( N ) )^(-1/theta[[g]][1]) - 1 )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- 1/theta[[g]][1]*( (1 - runif( N ) )^(-1/theta[[g]][1]) - 1 )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- 1/theta[[g]][1]*( (1 - runif( N ) )^(-1/1) - 1 )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <-             1/1*( (1 - runif( N ) )^(-1/theta[[g]][1]) - 1 )
    }
    if(family == "pstable")
    {
      w0     <- -log( runif(N) )
      theta0 <-  pi*( runif(N) - 1/2 )
      sigma0 <-  (cos(pi*theta[[g]][1]/4))^(2/theta[[g]][1])
      r0     <-  sin(theta[[g]][1]/2*(pi/2 + theta0))/( cos(theta[[g]][1]*pi/4)*cos(theta0) )^(2/theta[[g]][1])*
        ( cos(theta[[g]][1]*pi/4 + (theta[[g]][1]/2 - 1)*theta0)/w0 )^( (2 - theta[[g]][1])/theta[[g]][1] )
      W[, g] <- ( r0 - tan(pi*theta[[g]][1]/4) )*sigma0 + sigma0*tan(pi*theta[[g]][1]/4)
    }
    if(family == "ptstable") W[, g] <- 1/( 2^(1 - 2/theta[[g]][1])*rptstable( N, theta[[g]][1], Dim/2) )
    if(family == "rayleigh") W[, g] <- rweibull( N, shape = 2, scale = theta[[g]][1] )
    if(family == "weibull")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- rweibull( N, shape = theta[[g]][1], scale = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- rweibull( N, shape = theta[[g]][1], scale = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- rweibull( N, shape = 1            , scale = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <- rweibull( N, shape = theta[[g]][1], scale = 1             )
    }
    if(family == "constant") W[, g] <- rep(1, N)
    if( model != "restricted" )
    {
      if(family == "constant")
      {
        W[, g]    <- rep(1, N)
        T0[, , g] <- qnorm( ( runif(Q*N) + 1 )/2 )
        T1[, , g] <- sapply(1:N, function(r) T0[r, , g]*sqrt( W[r, g] ) )
      }else{
        T0[, , g] <- qnorm( ( runif(Q*N) + 1 )/2 )
        T1[, , g] <- sapply(1:N, function(r) T0[r, , g]*sqrt( W[r, g] ) )
      }
    }
    Mu     <- as.vector( mu[[g]] )
    Sigma  <- as.matrix ( sigma[[g]] )
    ############ computing f_y for g-th component of the restricted miture model                          ########
    if( model == "restricted" )
    {
      Lambda   <- as.vector( lambda[[g]] )
      Omega    <- Sigma + Lambda%*%t( Lambda )
      delta[g] <- 1 - mahalanobis( Lambda, rep(0, Dim), Omega )
      C0[g]    <- 2/sqrt( (2*pi)^( Dim + 1 )* abs( det( Sigma ) ) )
      dy[, g]  <- mahalanobis( Y, Mu, Omega )
      m[, g]   <- sapply(1:n, function(i) as.numeric(t( Lambda )%*%solve(Omega)%*%c(Y[i, ] - Mu)) )
      P1[g]    <- C0[g]*( 2*delta[g] )^( (0 + 1)/2 ) * gamma( (0 + 1)/2 )/2
      f_y[, g] <- weight[g]*P1[g]*colMeans( sapply(1:n, function(i){ W[, g]^( 1/2 - (Dim + 1)/2 )*exp( -0.5*dy[i, g]/W[, g] )*
                                                          (1 + pgamma( m[i, g]^2/(2*delta[g]*W[, g]), shape = 1/2,
                                                          rate = 1)*( sign( m[i, g] ) )^(0 + 1))
                                                         }
                                        )
                                )
      ############ end of computing f_y for g-th component of the restricted miture model                ########
    }else{
      Lambda       <- as.matrix( lambda[[g]] )
      ############ computing f_y for g-th component of the canonical or unrestricted miture model        ########
      C0[g]    <- 1/( (2*pi)^( Dim/2 )*sqrt( abs( det( Sigma ) ) ) )
      f_y[, g] <- weight[g]*sapply(1:n, function(i){ mean( exp( log( C0[g] ) + ( -Dim/2 )*log( W[, g] ) - 0.5*mahalanobis( t(
                                                     sapply(1:N, function(r) Y[i, ] - Mu - Lambda%*%T1[r, , g] ) ), rep(0, Dim ), Sigma )/W[, g] ) ) }                                  )
           ############ end of computing f_y for g-th component of the canonical or unrestricted miture model ########
    }
  }
  out <- apply( f_y, 1, sum )
  if(log == "TRUE") out <- log(out)
  return( out )
}
