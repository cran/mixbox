estep2 <- function(Y, G, weight, mu, sigma, lambda, family, skewness, param, theta, tick, h, N, PDF)
{
  Dim <- length( mu[[1]] )
  Q   <- length( lambda[[1]][1, ] )
  C0  <- rep( 0, G )
  P1  <- P2 <- P3 <- P4 <- rep( NA, N )
  Dim_theta <- length( theta[[1]] )
  theta1 <- vector( "list", G )
  theta2 <- vector( "list", G )
  n <- length( Y[, 1] )
  M <- floor( 0.7*N )
  K <- 4 + Dim_theta
  T2 <- array( NA, c( M, Dim_theta, Q ) )
  X1 <- X2 <- X3 <- X4 <- X5 <- matrix( NA, nrow = M, ncol = Dim )
  dy <- tau.hat <- my <- E_Zig <- E_H1_W0 <- f_y <- matrix( NA, nrow = n, ncol = G )
  E_deriv_theta <- array( NA, c( n, Dim_theta, G ) )
  Sigmainv <- list( )
  E_H1_T   <- array( 0, c( n, Q  , G ) )
  E_H1_TT  <- array( 0, c( n, Q*Q, G ) )
  if( family == "constant" )
  {
    H  <- array( NA, c( N, 4, G ) )
    T0 <- array( NA, c( N, 4, Q ) )
    T1 <- array( NA, c( N, 4, Q ) )
  }else{
    X0 <- rep(NA, N*K)
    H  <- array( NA, c( N, K, G ) )
    T0 <- array( NA, c( N, K, Q ) )
    T1 <- array( NA, c( N, K, Q ) )
  }
  for (g in 1:G)
  {
    if(family == "burriii")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- exp( -log( exp(- log( runif( N*K ) )/theta[[g]][1] ) - 1  )/theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- exp( -log( exp(- log( runif( N*K ) )/theta[[g]][1] ) - 1  )/theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- exp( -log( exp(- log( runif( N*K ) )/1             ) - 1  )/theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <- exp( -log( exp(- log( runif( N*K ) )/theta[[g]][1] ) - 1  )/1             )
    }
    if(family == "bs")
    {
      if( length(param) == 2 & sum( tick ) == 2 )
      {
        X0 <- theta[[g]][1]*rnorm(N*K)/2
        H[, , g] <- theta[[g]][2]*( 1 + 2*X0^2 + 2*X0*sqrt(1 + X0^2) )
      }
      if( length(param) == 1 & sum( tick ) == 2 )
      {
        X0 <- theta[[g]][1]*rnorm(N*K)/2
        H[, , g] <- theta[[1]][2]*( 1 + 2*X0^2 + 2*X0*sqrt(1 + X0^2) )
      }
      if( length(param) == 1 & tick[1] == 0 )
      {
        X0 <- 1*rnorm(N*K)/2
        H[, , g] <- theta[[1]][1]*( 1 + 2*X0^2 + 2*X0*sqrt(1 + X0^2) )
      }
      if( length(param) == 1 & tick[2] == 0 )
      {
        X0 <- theta[[g]][1]*rnorm(N*K)/2
        H[, , g] <- 1*( 1 + 2*X0^2 + 2*X0*sqrt(1 + X0^2) )
      }
    }
    if(family == "chisq")
    {
      H[, , g] <- rchisq(N*K, df = theta[[g]][1])
    }
    if(family == "f")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- rf(N*K, df1 = theta[[g]][1], df2 = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- rf(N*K, df1 = theta[[g]][1], df2 = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- rf(N*K, df1 = 1            , df2 = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <- rf(N*K, df1 = theta[[g]][1], df2 = 1             )
    }
    if(family == "gamma")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- rgamma(N*K, shape = theta[[g]][1], rate = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- rgamma(N*K, shape = theta[[g]][1], rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- rgamma(N*K, shape = 1            , rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <- rgamma(N*K, shape = theta[[g]][1], rate = 1             )
    }
    if(family == "igamma")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- 1/rgamma(N*K, shape = theta[[g]][1], rate = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- 1/rgamma(N*K, shape = theta[[g]][1], rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- 1/rgamma(N*K, shape = 1            , rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <- 1/rgamma(N*K, shape = theta[[g]][1], rate = 1             )
    }
    if(family == "gigaussian")
    {
      if( length(param) == 3 & sum( tick ) == 3 ) H[, , g] <- rgig( N*K, lambda = theta[[g]][1], chi = theta[[g]][2], psi = theta[[g]][3] )
      if( length(param) == 2 &     tick[3] == 0 ) H[, , g] <- rgig( N*K, lambda = theta[[g]][1], chi = theta[[g]][2], psi = theta[[g]][2] )
    }
    if(family == "igaussian")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- rigaussian( N*K, alpha = theta[[g]][1], beta = theta[[g]][2] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- rigaussian( N*K, alpha = 1            , beta = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <- rigaussian( N*K, alpha = theta[[g]][1], beta = 1             )
    }
    if(family == "lidley")
    {
      X0 <- rbinom( N*K, 1, theta[[g]][1]/(1 + theta[[g]][1]) )
      H[, , g] <- X0*rgamma( N*K, shape = 1, rate = theta[[g]][1] ) + (1 - X0)*rgamma( N*K, shape = 2, rate = theta[[g]][1] )
    }
    if(family == "loglog")
    {
      u <- runif( N*K )
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- theta[[g]][2]*( u/(1 - u) )^(1/theta[[g]][1] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- theta[[g]][1]*( u/(1 - u) )^(1/theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- theta[[g]][2]*( u/(1 - u) )^(1/1             )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <-             1*( u/(1 - u) )^(1/theta[[g]][1] )
    }
    if(family == "lognorm")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- rlnorm( N*K, meanlog = theta[[g]][1], sdlog = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- rlnorm( N*K, meanlog = theta[[g]][1], sdlog = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- rlnorm( N*K, meanlog = 0            , sdlog = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <- rlnorm( N*K, meanlog = theta[[g]][1], sdlog = 1             )
    }
    if(family == "lomax")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- 1/theta[[g]][2]*( (1 - runif( N*K ) )^(-1/theta[[g]][1]) - 1 )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- 1/theta[[g]][1]*( (1 - runif( N*K ) )^(-1/theta[[g]][1]) - 1 )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- 1/theta[[g]][1]*( (1 - runif( N*K ) )^(-1/1) - 1 )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <-             1/1*( (1 - runif( N*K ) )^(-1/theta[[g]][1]) - 1 )
    }
    if(family == "pstable")
    {
      w0     <- -log( runif(N*K) )
      theta0 <-  pi*( runif(N*K) - 1/2 )
      sigma0 <-  (cos(pi*theta[[g]][1]/4))^(2/theta[[g]][1])
      r0     <-  sin(theta[[g]][1]/2*(pi/2 + theta0))/( cos(theta[[g]][1]*pi/4)*cos(theta0) )^(2/theta[[g]][1])*
        ( cos(theta[[g]][1]*pi/4 + (theta[[g]][1]/2 - 1)*theta0)/w0 )^( (2 - theta[[g]][1])/theta[[g]][1] )
      H[, , g] <-  ( r0 - tan(pi*theta[[g]][1]/4) )*sigma0 + sigma0*tan(pi*theta[[g]][1]/4)
    }
    if(family == "ptstable")
    {
      H[, , g] <- 1/( 2^(1 - 2/theta[[g]][1])*rptstable(N*K, theta[[g]][1], Dim/2) )
    }
    if(family == "rayleigh")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- rweibull(N*K, shape = theta[[g]][1], scale = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- rweibull(N*K, shape = theta[[g]][1], scale = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- rweibull(N*K, shape = 1            , scale = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <- rweibull(N*K, shape = theta[[g]][1], scale = 1             )
    }
    if(family == "weibull")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) H[, , g] <- rweibull(N*K, shape = theta[[g]][1], scale = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) H[, , g] <- rweibull(N*K, shape = theta[[g]][1], scale = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) H[, , g] <- rweibull(N*K, shape = 1            , scale = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) H[, , g] <- rweibull(N*K, shape = theta[[g]][1], scale = 1             )
    }
    if(family == "constant")
    {
      H[, , g]    <- runif(N*4, 0.9999, 1.0001)
      T0[, 1:4, 1:Q] <- qnorm( ( runif(Q*N*4) + 1 )/2 )
      T1[, 1:4, 1:Q] <- sweep( T0[, 1:4, 1:Q], c(1, 2), sqrt( H[, 1:4, g] ), "*" )
    }else{
      T0[, 1:4, 1:Q] <- qnorm( ( runif(Q*N*4) + 1 )/2 )
      T1[, 1:4, 1:Q] <- sweep( T0[, 1:4, 1:Q], c(1, 2), sqrt( H[, 1:4, g] ), "*" )
    }
    ############ computing  f_y_g ########
    Mu     <- as.vector( mu[[g]] )
    Sigma  <- as.matrix( sigma[[g]] )
    Lambda <- as.matrix( lambda[[g]] )
    C0[g] <- 1/( (2*pi)^( Dim/2 )*sqrt( abs( det( Sigma ) ) ) )
    for(i in 1:n)
    {
      ############ computing f(y)                      ########
      index <- sample( 1:N, M )
      T2[, 1, 1:Q] <- T1[index, 1, 1:Q]
      X1 <- t( sapply(1:M, function(r) Y[i, ] - Mu - Lambda%*%T2[r, 1, 1:Q] ) )
      P1 <- exp( log( C0[g] ) + ( -Dim/2     )*log( H[index, 1, g] ) -0.5*mahalanobis( X1, rep(0, Dim ), Sigma )/H[index, 1, g] )
      f_y[i, g] <- mean( P1 )
      ############ end of computing f(y)               ########
      ############ computing E(W^(-1) | Y)            ########
      index <- sample( 1:N, M )
      T2[, 1, 1:Q] <- T1[index, 2, 1:Q]
      X2 <- t( sapply(1:M, function(r) Y[i, ] - Mu - Lambda%*%T2[r, 1, 1:Q] ) )
      P2 <- exp( log( C0[g] ) + ( -Dim/2 - 1 )*log( H[index, 2, g] ) -0.5*mahalanobis( X2, rep(0, Dim ), Sigma )/H[index, 2, g] )
      E_H1_W0[i, g] <- mean( P2 )/f_y[i, g]
      ############ end of computing E(ZW^(-1) | Y)     ########
      ############ computing E(W^(-1)T | Y)           ########
      index <- sample( 1:N, M )
      T2[, 1, 1:Q] <- T1[index, 3, 1:Q]
      X3 <- t( sapply(1:M, function(r) Y[i, ] - Mu - Lambda%*%T2[r, 1, 1:Q] ) )
      P3 <- exp( log( C0[g] ) + ( -Dim/2 - 1 )*log( H[index, 3, g] ) -0.5*mahalanobis( X3, rep(0, Dim ), Sigma )/H[index, 3, g] )
      E_H1_T[i, , g] <- apply( sweep( T2[, 1, 1:Q], 1, P3, "*" ), 2, mean, na.rm = TRUE )/f_y[i, g]
      ############ end of computing E(ZW^(-1)T | Y)    ########
      ############ computing E(W^(-1)T^tT | Y)        ########
      index <- sample( 1:N, M )
      T2[, 1, 1:Q] <- T1[index, 4, 1:Q]
      X4 <- t( sapply(1:M, function(j) Y[i, ] - Mu - Lambda%*%T2[j, 1, 1:Q] ) )
      P4 <- exp( log( C0[g] ) + ( -Dim/2 - 1 )*log( H[index, 4, g] ) -0.5*mahalanobis( X4, rep(0, Dim ), Sigma )/H[index, 4, g] )
      E_H1_TT[i, , g] <- matrix( apply( sapply(1:M, function(r) T2[r, 1, 1:Q]%o%T2[r, 1, 1:Q]*P4[r] ), 1, mean, na.rm = TRUE ),
                                 nrow = Q, ncol = Q )/f_y[i, g]
      ############ end of computing E(W^(-1)T^tT | Y) ########
    }
  }
  ############ computing E(Z | Y)        ########
  E_Zig <- sweep( f_y, 2, weight, "*" )
  E_Zig <- E_Zig/rowSums( E_Zig )
  ############ end of computing E(Z | Y) ########
  for( g in 1:G )
  {
    Mu     <- as.vector( mu[[g]] )
    Sigma  <- as.matrix( sigma[[g]] )
    Lambda <- as.matrix( lambda[[g]] )
    E_H1_W0[, g] <- E_H1_W0[, g]*E_Zig[, g]
    E_H1_T[, , g] <- E_H1_T[, , g]*E_Zig[, g]
    E_H1_TT[, , g] <- E_H1_TT[, , g]*E_Zig[, g]
    ############ computing derivative of f_W(w) w.r.t theta ########
    if( any( family == c( "pstable", "ptstable" ) ) )
    {
      index <- sample( 1:N, M )
      T0[, 5, 1:Q] <- qnorm( ( runif(Q*N) + 1 )/2 )
      T1[, 5, 1:Q] <- T0[, 5, 1:Q]*sqrt( H[, 5, g] )
    if( family == "pstable")
    {
      forward_diff  <- dpstable( H[index, 5, g], theta[[g]] + h )
      backward_diff <- dpstable( H[index, 5, g], theta[[g]] - h )
      derive_theta  <- ( forward_diff - backward_diff )/( 2*h*dpstable( H[index, 5, g], theta[[g]] ) )
    }else{
      forward_diff  <- dptstable( H[index, 5, g], theta[[g]] + h )
      backward_diff <- dptstable( H[index, 5, g], theta[[g]] - h )
      derive_theta  <- ( forward_diff - backward_diff )/( 2*h*dptstable( H[index, 5, g], theta[[g]] ) )
    }
      for(i in 1:n)
      {
        index <- sample( 1:N, M )
        T2[, 1, 1:Q] <- T1[index, 5, 1:Q]
        X5 <- t( sapply(1:M, function(r) Y[i, ] - Mu - Lambda%*%T2[r, 1, 1:Q] ) )
        P5 <- exp( log( C0[g] ) -( Dim/2 )*log( H[index, 5, g] ) - 0.5*mahalanobis( X5, rep(0, Q ), Sigma )/H[index, 5, g] )
        E_deriv_theta[i, 1, g] <- mean( derive_theta*P5, na.rm = TRUE)/f_y[i, g]
      }
    }
    if( all( family != c( "pstable", "ptstable", "constant" ) ) )
    {
      T0[, 5: (4 + Dim_theta), 1:Q] <- qnorm( ( runif(Q*N*Dim_theta) + 1 )/2 )
      T1[, 5: (4 + Dim_theta), 1:Q] <- sweep( T0[, 5: (4 + Dim_theta), 1:Q], c(1, 2), sqrt( H[, 5: (4 + Dim_theta), g] ), "*" )
      first_deriv <- tryCatch( sapply( 1:Dim_theta, function(r) D( bquote( .( PDF ) ), param[r] ) ), error = function(e)( "fail" )  )
      if( any( first_deriv == "fail" ) )
      {
        pdf1 <- function(x, param){ }
        body( pdf1 ) <- bquote( .( PDF ) )
        theta2[[g]] <- theta[[g]]
        theta1[[g]] <- theta[[g]]
        for(j in 1:Dim_theta)
        {
          theta2[[g]][j] <- theta[[g]][j] + h
          theta1[[g]][j] <- theta[[g]][j] - h
          for(k in 1:Dim_theta) assign( param[k], theta2[[g]][k] )
          forward_diff <- pdf1( H[index, 4 + j, g], theta2[[g]] )
           for(k in 1:Dim_theta) assign( param[k], theta1[[g]][k] )
            backward_diff <- pdf1( H[index, 4 + j, g], theta1[[g]] )
            for(k in 1:Dim_theta) assign( param[k], theta[[g]][k] )
            derive_theta <- ( forward_diff - backward_diff )/( 2*h*pdf1( H[index, 4 + j, g], theta[[g]] ) )
          for(i in 1:n)
          {
            index <- sample( 1:N, M )
            T2[, (1:Dim_theta), 1:Q] <- T1[index, 5: (4 + Dim_theta), 1:Q]
            X5 <- t( sapply(1:M, function(r) Y[i, ] - Mu - Lambda%*%T2[r, j, 1:Q] ) )
            P5 <- exp( log( C0[g] ) -( Dim/2 )*log( H[index, 4 + j, g] ) - 0.5*mahalanobis( X5, rep(0, Q ), Sigma )/H[index, 4 + j, g] )
            E_deriv_theta[i, j, g] <- mean( derive_theta*P5, na.rm = TRUE)/f_y[i, g]
          }
          theta1[[g]] <- theta[[g]]
          theta2[[g]] <- theta[[g]]
        }
      }else{
        pdf2 <- function(x, y, i){ }
        for(k in 1:Dim_theta) assign( param[k], theta[[g]][k] )
        for(j in 1:Dim_theta)
        {
          body(pdf2) <- bquote(	.( first_deriv[[j]] )/.( PDF )*.( quote( exp( log( C0[g] ) -( Dim/2 )*log( x ) -
                                0.5*mahalanobis( y, rep(0, Q ), Sigma )/x ) ) )	)
          for(i in 1:n)
          {
            index <- sample( 1:N, M )
            T2[, (1:Dim_theta), 1:Q] <- T1[index, 5: (4 + Dim_theta), 1:Q]
            X5 <- t( sapply(1:M, function(r) Y[i, ] - Mu - Lambda%*%T2[r, j, 1:Q] ) )
            E_deriv_theta[i, j, g] <- mean( pdf2( H[index, 4 + j, g], X5, i), na.rm = TRUE )/f_y[i, g]
          }
        }
      }
    }
    E_deriv_theta[, , g] <- E_deriv_theta[, , g]*E_Zig[, g]
  }
  ########### end of computing derivative of f_W(w) w.r.t theta ########
  for(g in 1:G)
  {
    for(i in 1:n)
    {
      if( is.nan( E_Zig[i, g] ) ) E_Zig[i, g]  <- runif(1, 0, .Machine$double.xmin)
      if( is.nan( E_H1_W0[i, g] ) | is.infinite( E_H1_W0[i, g] ) )
      {
        E_H1_W0[i, g] <- runif(1, 0, .Machine$double.xmin)
      }
      for(j in 1:Dim_theta)
      {
        if( is.nan( E_deriv_theta[i, j, g] ) ) E_deriv_theta[i, j, g] <- 0
        if( is.infinite( E_deriv_theta[i, j, g] ) ) E_deriv_theta[i, j, g] <- 0
      }
      for(k in 1:Q)
      {
        if( is.nan( E_H1_T[i, k, g] ) | is.infinite( E_H1_T[i, k, g] ) )
        {
          E_H1_T[i, k, g] <- runif(1, 0, .Machine$double.xmin)
        }
      }
      A <- matrix( E_H1_TT[i, , g], nrow = Q, ncol = Q )
      for(r in 1:Q)
      {
        for(s in 1:Q)
        {
          if( is.nan( A[r, s] ) | is.infinite( A[r, s] ) )
          {
            A[r, s] <- runif(1, 0, .Machine$double.xmin)
          }
        }
      }
      E_H1_TT[i, , g] <- ( A + t(A) )/2
    }
  }
  return( list(f_y = f_y, E_Zig = E_Zig, E_H1_W0 = E_H1_W0, E_deriv_theta = E_deriv_theta, E_H1_TT = E_H1_TT,
               E_H1_T = E_H1_T) )
}
