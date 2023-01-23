sefm <- function(Y, G, weight, model = "restricted", mu, sigma, lambda, family = "constant", skewness = "FALSE", param = NULL, theta = NULL, tick = rep(1, 2), h = 0.001, N = 3000, level = 0.05, PDF = NULL)
{
  if( sum( weight ) > 1 || any( weight < 0 ) ) stop( "Elements of mixing proportion must be positive and must sum to 1." )
  if( length( weight ) != G ) stop( "Length of mixing proportions and number of components G must be equal." )
  if( length( mu     ) != G ) stop( "Length of the location parameter and number of components G must be equal." )
  if( length( sigma  ) != G ) stop( "Length of the dispersion matrix and number of components G must be equal." )
  if( length( lambda ) != G ) stop( "Length of the skewness parameter and number of components G must be equal." )
  if(family != "constant" & family != "bs" & family != "burriii" & family != "chisq" & family != "exp" & family != "f" & family != "gamma" & family != "gigaussian"  & family != "igamma" & family != "lidley" &
     family != "loglog" & family != "lognorm" & family != "lomax" & family != "pstable" & family != "ptstable" & family != "rayleigh" & family != "weibull" )
  { stop( "Mixing distribution not implemented or misspelled. Please check the manual for guidelines." ) }
  if( skewness != TRUE & skewness != FALSE ) stop( "Skewness must be a logical statement either TRUE or FALSE." )
  if( family == "constant" )
  {
    param <- c( "nu" )
    theta <- rep( list(200), G )
    tick  <- rep(1, 1)
    h     <- 0.001
    PDF   <- quote( (nu/2)^(nu/2)*x^(nu/2 - 1)/gamma(nu/2)*exp(-nu*x/2) )
  }
  if( family != "constant" & family != "pstable" & family != "ptstable" & is.null( PDF ) ) stop( "Expression for the density function of mixing distribution must be given." )
  if( family != "constant" & is.null( param ) ) stop( "Name of the mixing distribution parameters must be given." )
  if( family != "constant" & is.null( theta ) ) stop( "ML estimator of the mixing distribution parameters must be given." )
  if( family != "constant" & is.null( tick ) ) stop( " vector tick must be given." )
  if( family != "constant" & length( theta ) != G ) stop( "Length of the ML estimators of mixing distribution and number of components G must be equal." )
  if( length( param ) != length( theta[[1]] ) ) stop( "Length of the parameter vector of mixing distribution and associated MLEs must be equal." )
  if( any( tick < 0 ) || any( tick > 1 ) || sum( tick ) == 0 || ( sum( tick )%%1 != 0 & exp( prod( tick ) ) != exp(1) ) ) stop( "Elements of vector tick are either 0 or 1." )
  if( length(tick) < length(param) ) stop( "Length of vector tick cannot exceed the length of param." )
  Dim <- length( Y[1, ] )
  if( all( model != c("canonical", "restricted", "unrestricted") ) ) stop( "model must be canonical, restricted, or unrestricted." )
  if(model == "restricted") Q <- length( lambda[[1]] )
  if(model == "canonical" | model == "unrestricted" )  Q <- dim( lambda[[1]] )[2]
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
  if( sum( s1 ) == 0 ) stop( "Skewnesss vector must be non-zero." )
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
              if( all( as.matrix(lambda[[g]])[i, ] == 0 ) ) s2[g] <- 0
            }
              for(j in 1:Q)
              {
                if( all( as.matrix(lambda[[g]])[ , j] == 0 ) ) s3[g] <- 0
              }
          }
        }else{
          lambda <- vector("list", G)
            for(g in 1:G)
            {
              lambda[[g]] <- matrix( 0, nrow = Dim, ncol = Q )
            }
        }
      if( sum( s2 ) == 0 | sum( s3 ) == 0  ) stop( "Skewnesss vector must be non-zero." )
    }
      if(model == "restricted")
      {
        out_ofim1 <- ofim1(Y = Y, G = G, weight = weight, mu = mu, sigma = sigma, lambda = lambda, family = family, skewness = skewness, param = param, theta =
                       theta, tick = tick, h = h, N = N, level = level, PDF = PDF)
        out <- configuration1(Y = Y, G = G, weight = weight, mu = mu, sigma = sigma, lambda = lambda, family = family, skewness = skewness, param = param, theta =
                           theta, ofim1 = out_ofim1, level = level)
      }else{
        out_ofim2 <- ofim2(Y = Y, G = G, weight = weight, model = model, mu = mu, sigma = sigma, lambda = lambda, family = family, skewness = skewness, param = param, theta =
                     theta, tick = tick, h = h, N = N, level = level, PDF = PDF)
        out <- configuration2(Y = Y, G = G, weight = weight, model = model, mu = mu, sigma = sigma, lambda = lambda, family = family, skewness = skewness, param = param, theta =
                         theta, ofim2 = out_ofim2, level = level)
      }
  out
}
