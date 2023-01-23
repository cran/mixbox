estep1 <- function(Y, G, weight, mu, sigma, lambda, family, skewness, param, theta, tick, h, N, PDF)
{
  Dim <- length( mu[[1]] )
  P1 <- C0 <- rep( 0, G )
  Dim_theta <- length( theta[[1]] )
  theta2  <- theta1 <- vector( "list", G )
  n <- length( Y[, 1] )
  M <- floor( 0.7*N )
  delta <- rep( 0, G )
  dy <- tau.hat <- my <- E_Zig <- E_H1_W0 <- f_y <- matrix( NA, nrow = n, ncol = G )
  E_deriv_theta <- array( 0, c( n, Dim_theta, G ) )
  Sigmainv <- list( )
  E_H1_Wg  <- e_1j_g <- array( 0, c( n, 2, G ) )
  if( family == "constant" )
  {
    H <- array( runif(N*4*G, 0.9999, 1.0001), c(N, 4, G) )
  }else{
    H <- array( NA, c( N, 5, G ) )
  }
  for (g in 1:G)
  {
    if(family == "burriii")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- exp( -log( exp(- log( runif( N*5 ) )/theta[[g]][1] ) - 1  )/theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- exp( -log( exp(- log( runif( N*5 ) )/theta[[g]][1] ) - 1  )/theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- exp( -log( exp(- log( runif( N*5 ) )/1             ) - 1  )/theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <- exp( -log( exp(- log( runif( N*5 ) )/theta[[g]][1] ) - 1  )/1             )
    }
    if(family == "bs")
    {
      if( length(param) == 2 & sum( tick ) == 2 )
      {
        z0 <- theta[[g]][1]*rnorm(N*5)/2
        H[, , g] <- theta[[g]][2]*( 1 + 2*z0^2 + 2*z0*sqrt(1 + z0^2) )
      }
      if( length(param) == 1 & sum( tick ) == 2 )
      {
        z0 <- theta[[g]][1]*rnorm(N*5)/2
        H[, , g] <- theta[[1]][2]*( 1 + 2*z0^2 + 2*z0*sqrt(1 + z0^2) )
      }
      if( length(param) == 1 & tick[1] == 0 )
      {
        z0 <- 1*rnorm(N*5)/2
        H[, , g] <- theta[[1]][1]*( 1 + 2*z0^2 + 2*z0*sqrt(1 + z0^2) )
      }
      if( length(param) == 1 & tick[2] == 0 )
      {
        z0 <- theta[[g]][1]*rnorm(N*5)/2
        H[, , g] <- 1*( 1 + 2*z0^2 + 2*z0*sqrt(1 + z0^2) )
      }
    }
    if(family == "chisq")
    {
      H[, , g] <- rchisq(N*5, df = theta[[g]][1])
    }
    if(family == "f")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- rf(N*5, df1 = theta[[g]][1], df2 = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- rf(N*5, df1 = theta[[g]][1], df2 = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- rf(N*5, df1 = 1            , df2 = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <- rf(N*5, df1 = theta[[g]][1], df2 = 1             )
    }
    if(family == "gamma")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- rgamma(N*5, shape = theta[[g]][1], rate = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- rgamma(N*5, shape = theta[[g]][1], rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- rgamma(N*5, shape = 1            , rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <- rgamma(N*5, shape = theta[[g]][1], rate = 1             )
    }
    if(family == "igamma")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- 1/rgamma(N*5, shape = theta[[g]][1], rate = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- 1/rgamma(N*5, shape = theta[[g]][1], rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- 1/rgamma(N*5, shape = 1            , rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <- 1/rgamma(N*5, shape = theta[[g]][1], rate = 1             )
    }
    if(family == "gigaussian")
    {
      if( length(param) == 3 & sum( tick ) == 3 ) H[, , g] <- rgig( N*5, lambda = theta[[g]][1], chi = theta[[g]][2], psi = theta[[g]][3] )
      if( length(param) == 2 &     tick[3] == 0 ) H[, , g] <- rgig( N*5, lambda = theta[[g]][1], chi = theta[[g]][2], psi = theta[[g]][2] )
    }
    if(family == "igaussian")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- rigaussian( N*5, alpha = theta[[g]][1], beta = theta[[g]][2] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- rigaussian( N*5, alpha = 1            , beta = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <- rigaussian( N*5, alpha = theta[[g]][1], beta = 1             )
    }
    if(family == "lidley")
    {
      x <- rep( NA, N*5 )
      x <- rbinom( N*5, 1, theta[[g]][1]/(1 + theta[[g]][1]) )
      H[, , g] <- x*rgamma( N*5, shape = 1, rate = theta[[g]][1] ) + (1 - x)*rgamma( N*5, shape = 2, rate = theta[[g]][1] )
    }
    if(family == "loglog")
    {
      u <- runif( N*5 )
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- theta[[g]][2]*( u/(1 - u) )^(1/theta[[g]][1] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- theta[[g]][1]*( u/(1 - u) )^(1/theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- theta[[g]][2]*( u/(1 - u) )^(1/1             )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <-             1*( u/(1 - u) )^(1/theta[[g]][1] )
    }
    if(family == "lognorm")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- rlnorm( N*5, meanlog = theta[[g]][1], sdlog = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- rlnorm( N*5, meanlog = theta[[g]][1], sdlog = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- rlnorm( N*5, meanlog = 0            , sdlog = theta[[g]][1] )
      if( length(param) == 1 &    tick[2] == 0 ) H[, , g] <- rlnorm( N*5, meanlog = theta[[g]][1], sdlog = 1             )
    }
    if(family == "lomax")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- 1/theta[[g]][2]*( (1 - runif( N*5 ) )^(-1/theta[[g]][1]) - 1 )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- 1/theta[[g]][1]*( (1 - runif( N*5 ) )^(-1/theta[[g]][1]) - 1 )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- 1/theta[[g]][1]*( (1 - runif( N*5 ) )^(-1/1) - 1 )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <-             1/1*( (1 - runif( N*5 ) )^(-1/theta[[g]][1]) - 1 )
    }
    if(family == "pstable")
    {
      w0       <- -log( runif(N*5) )
      theta0   <-  pi*( runif(N*5) - 1/2 )
      sigma0   <-  (cos(pi*theta[[g]][1]/4))^(2/theta[[g]][1])
      r0       <-  sin(theta[[g]][1]/2*(pi/2 + theta0))/( cos(theta[[g]][1]*pi/4)*cos(theta0) )^(2/theta[[g]][1])*
        ( cos(theta[[g]][1]*pi/4 + (theta[[g]][1]/2 - 1)*theta0)/w0 )^( (2 - theta[[g]][1])/theta[[g]][1] )
      H[, , g] <-  ( r0 - tan(pi*theta[[g]][1]/4) )*sigma0 + sigma0*tan(pi*theta[[g]][1]/4)
        #rstable(N*5,theta[[g]][1]/2,1,cos(pi*theta[[g]][1]/4)^(2/theta[[g]][1]),0,1)
    }
    if(family == "ptstable")
    {
      H[, , g] <- 1/( 2^(1 - 2/theta[[g]][1])*rptstable(N*5, theta[[g]][1], Dim/2) )
    }
    if(family == "rayleigh")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- rweibull(N*5, shape = theta[[g]][1], scale = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- rweibull(N*5, shape = theta[[g]][1], scale = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- rweibull(N*5, shape = 1            , scale = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <- rweibull(N*5, shape = theta[[g]][1], scale = 1             )
    }
    if(family == "weibull")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- rweibull(N*5, shape = theta[[g]][1], scale = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- rweibull(N*5, shape = theta[[g]][1], scale = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- rweibull(N*5, shape = 1            , scale = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <- rweibull(N*5, shape = theta[[g]][1], scale = 1             )
    }
    #    }
    ############ computing  f_y_g ########
    Mu     <- as.vector( mu[[g]] )
    Sigma  <- as.matrix ( sigma[[g]] )
    Lambda <- as.vector( lambda[[g]] )
    Omega  <- Sigma + Lambda%*%t( Lambda )
    delta[g] <- 1 - mahalanobis( Lambda, rep(0, Dim), Omega )
    Sigmainv[[g]] <- solve( Sigma )
    C0[g] <- 2/sqrt( (2*pi)^( Dim + 1 )* abs( det( Sigma ) ) )
    dy[, g] <- mahalanobis( Y, Mu, Omega )
    my[, g] <- sapply(1:n, function(i) as.numeric(t( Lambda )%*%solve(Omega)%*%c(Y[i, ] - Mu)) )
    r  <- 0
    P1[g] <- C0[g]*( 2*delta[g] )^( (0 + 1)/2 ) * gamma( (0 + 1)/2 )/2
    X <- sample( H[, 1, g], M )
    P2 <- colMeans( sapply(1:n, function(i){ X^( (0 + 1)/2 - (Dim + 1)/2 + r)*exp( -0.5*dy[i, g]/X )*(1 +
                                             pgamma( my[i, g]^2/(2*delta[g]*X), shape = (0 + 1)/2, rate = 1)*(-1)^(0)*( sign(my[i, g]) )^(0 + 1))} ) )
    f_y[, g] <- P1[g]*P2
    ############ end of computing  f(y)   ########
    ############ computing E(ZW^(-1) | Y) ########
    r <- -1
    X <- sample( H[, 2, g], M )
    P2 <- colMeans( sapply(1:n, function(i){ X^( (0 + 1)/2 - (Dim + 1)/2 + r)*exp( -0.5*dy[i, g]/X )*(1 +
                                             pgamma( my[i, g]^2/(2*delta[g]*X), shape = (0 + 1)/2, rate = 1)*(-1)^(0)*( sign(my[i, g]) )^(0 + 1))} ) )
    E_H1_W0[, g] <- P1[g]*P2/f_y[, g]
    ############ end of computing  E(ZW^(-1) | Y) ########
    ############ computing E(ZW^(-1)T | Y)        ########
    for( j in 1:2 )
    {
      X <- sample( H[, (j + 2), g], M )
      P3 <- C0[g]*( 2*delta[g] )^( (j + 1)/2 ) * gamma( (j + 1)/2 )/2
      P2 <- colMeans( sapply(1:n, function(i){ X^( (j + 1)/2 - (Dim + 1)/2 + r)*exp( -0.5*dy[i, g]/X )*(1 +
                                               pgamma( my[i, g]^2/(2*delta[g]*X), shape = (j + 1)/2, rate = 1)*(-1)^(j)*( sign(my[i, g]) )^(j + 1)) } ) )
      e_1j_g[, j,  g] <- P3*P2/f_y[, g]
    }
    E_H1_Wg[, 1, g] <- (my[, g]*E_H1_W0[, g] + e_1j_g[, 1, g] )
    E_H1_Wg[, 2, g] <- (my[, g]^2*E_H1_W0[, g] + 2*my[, g]*e_1j_g[, 1, g] + e_1j_g[, 2, g] )
  }
  ############ end of computing E(ZW^(-1)T | Y) ########
  ############ computing E(Z | Y)               ########
  E_Zig <- sweep( f_y, 2, weight, "*" )
  E_Zig <- E_Zig/rowSums( E_Zig )
  ############ end of computing E(Z | Y) ########
  for( g in 1:G )
  {
    E_H1_W0[,    g] <- E_H1_W0[,   g]*E_Zig[, g]
    E_H1_Wg[, ,  g] <- E_H1_Wg[, , g]*E_Zig[, g]
    ############ computing derivative of f_W(w) w.r.t theta ########
    if( all( family != c( "pstable", "ptstable", "constant" ) ) )
    {
      first_deriv <- tryCatch( sapply( 1:Dim_theta, function(k) D( bquote( .( PDF ) ), param[k] ) ), error = function(e)( "fail" )  )
      pdf1 <- function(x, param){ }
      body( pdf1 ) <- bquote( .( PDF ) )
      X <- sample( H[, 5, g], M )
      if( any(first_deriv == "fail") )
      {
        theta2[[g]] <- theta[[g]]
        theta1[[g]] <- theta[[g]]
        for(j in 1:Dim_theta)
        {
          theta2[[g]][j] <- theta[[g]][j] + h
          theta1[[g]][j] <- theta[[g]][j] - h
            for(k in 1:Dim_theta) assign( param[k], theta2[[g]][k] )
            forward_diff <- pdf1( X, theta2[[g]] )
            for(k in 1:Dim_theta) assign( param[k], theta1[[g]][k] )
            backward_diff <- pdf1( X, theta1[[g]] )
            for(k in 1:Dim_theta) assign( param[k], theta[[g]][k] )
            derive_theta <- ( forward_diff - backward_diff )/( 2*h*pdf1( X, theta[[g]] ) )
            for(i in 1:n)
            {
            E_deriv_theta[i, j, g] <- 2*P1[g]*mean( derive_theta*X^(-Dim/2)*exp(-dy[i, g]/(2*X))*
                                                    pnorm( my[i, g] /sqrt( X*delta[g] ) ), na.rm = TRUE)/f_y[i, g]
          }
          theta1[[g]] <- theta[[g]]
          theta2[[g]] <- theta[[g]]
        }
      }else{
        pdf2 <- function(x, i){ }
        for(k in 1:Dim_theta) assign( param[k], theta[[g]][k] )
        for(j in 1:Dim_theta)
        {
          body(pdf2) <- bquote( .( first_deriv[[j]] )/.( PDF ) *.( quote( x^(-Dim/2) ) ) *.( quote( exp( -dy[i, g]/(2*x) )) ) *
                                .( quote( pnorm( ( my[i, g] )/sqrt(x*delta[g]) ) ) ) )
          for(i in 1:n)
          {
            E_deriv_theta[i, j, g] <- 2*P1[g]*mean( pdf2( X, i), na.rm = TRUE )/f_y[i, g]
          }
        }
      }
    }
    if( family == "pstable")
    {
        X <- sample( H[, 5, g], M )
        forward_diff  <- dpstable( X, theta[[g]] + h )
        backward_diff <- dpstable( X, theta[[g]] - h )
        derive_theta  <- ( forward_diff - backward_diff )/( 2*h*dpstable( X, theta[[g]] ) )
        for(i in 1:n)
        {
          E_deriv_theta[i, , g] <- 2*P1[g]*mean( derive_theta*X^(-Dim/2)*exp( -dy[i, g]/(2*X) )*
                                                 pnorm( my[i, g] /sqrt( X*delta[g] ) ), na.rm = TRUE)/f_y[i, g]
        }
    }
    if( family == "ptstable")
    {
          X <- sample( H[, 5, g], M )
          forward_diff  <- dptstable( X, theta[[g]] + h )
          backward_diff <- dptstable( X, theta[[g]] - h )
          derive_theta  <- ( forward_diff - backward_diff )/( 2*h*dptstable( X, theta[[g]] ) )
          for(i in 1:n)
          {
          E_deriv_theta[i, , g] <- 2*P1[g]*mean( derive_theta*X^(-Dim/2)*exp(-dy[i, g]/(2*X))*
                                                 pnorm( my[i, g] /sqrt( X*delta[g] ) ), na.rm = TRUE)/f_y[i, g]
          }
    }
    E_deriv_theta[, , g] <- E_deriv_theta[, , g]*E_Zig[, g]
  }
  ############ end of computing derivative of f_W(w) w.r.t theta ########
  for(g in 1:G)
  {
    for(i in 1:n)
    {
      for(j in 1:Dim_theta)
      {
        if( is.nan(E_deriv_theta[i, j, g]) == TRUE ) E_deriv_theta[i, j, g] <- 0
        if( is.infinite(E_deriv_theta[i, j, g]) == TRUE) E_deriv_theta[i, j, g] <- 0
      }
      if( is.nan(E_Zig[i, g]) || is.nan(E_H1_W0[i, g]) || is.nan(E_H1_Wg[i, 1, g]) || is.nan(E_H1_Wg[i, 2, g]) )
      {
        E_Zig[i, g] <- runif(1, 0, .Machine$double.xmin);     E_H1_W0[i, g] <- runif(1, 0, .Machine$double.xmin);
        E_H1_Wg[i, 1, g]<- runif(1, 0, .Machine$double.xmin); E_H1_Wg[i, 2, g]<- runif(1, 0, .Machine$double.xmin);
      }
    }
  }
  return( list(f_y = f_y, E_Zig = E_Zig, E_H1_W0 = E_H1_W0, E_deriv_theta = E_deriv_theta, E_H1_Wg = E_H1_Wg, my = my, dy = dy, delta = delta,
               Sigmainv = Sigmainv, Dim = Dim) )
}

