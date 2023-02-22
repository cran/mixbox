estep2 <- function(Y, G, weight, mu, sigma, lambda, family, skewness, param, theta, tick, h, N, PDF)
{
  Dim <- length( mu[[1]] )
  Q   <- length( lambda[[1]][1, ] )
  C0  <- rep( 0, G )
  P1  <- P2 <- P3 <- P4 <- rep( NA, N )
  theta1 <- vector( "list", G )
  theta2 <- vector( "list", G )
  n <- length( Y[, 1] )
  M <- N #floor( 0.9*N )
  X1 <- X2 <- X3 <- X4 <- X5 <- matrix( NA, nrow = M, ncol = Dim )
  X6 <- matrix( NA, nrow = M, ncol = Q )
  dy <- tau_ig <- my <- tau_ig <- E_W1 <- f_y <- matrix( NA, nrow = n, ncol = G )
  if( family == "constant" )
  {
    Dim_theta <- 1
    W  <- matrix( NA, nrow = N*4, ncol = G )
    T0 <- array( NA, c( N*4, Q, G ) )
    T1 <- array( NA, c( N*4, Q, G ) )
  }else{
    Dim_theta <- length( theta[[1]] )
    K <- 4 + Dim_theta
    X0 <- rep(NA, N*K)
    W  <- matrix( NA, nrow = N*K, ncol = G )
    T0 <- array( NA, c( N*K, Q, G ) )
    T1 <- array( NA, c( N*K, Q, G ) )
  }
  E_deriv_theta <- array( NA, c( n, Dim_theta, G ) )
  Sigmainv <- list( )
  E_W1_T   <- array( 0, c( n, Q  , G ) )
  E_W1_TT  <- array( 0, c( n, Q*Q, G ) )
  for (g in 1:G)
  {
    if(family == "burriii")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- exp( -log( exp(- log( runif( N*K ) )/theta[[g]][1] ) - 1  )/theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- exp( -log( exp(- log( runif( N*K ) )/theta[[g]][1] ) - 1  )/theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- exp( -log( exp(- log( runif( N*K ) )/1             ) - 1  )/theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <- exp( -log( exp(- log( runif( N*K ) )/theta[[g]][1] ) - 1  )/1             )
    }
    if(family == "bs")
    {
      if( length(param) == 2 & sum( tick ) == 2 )
      {
        X0 <- theta[[g]][1]*rnorm(N*K)/2
        W[, g] <- theta[[g]][2]*( 1 + 2*X0^2 + 2*X0*sqrt(1 + X0^2) )
      }
      if( length(param) == 1 & sum( tick ) == 2 )
      {
        X0 <- theta[[g]][1]*rnorm(N*K)/2
        W[, g] <- theta[[1]][2]*( 1 + 2*X0^2 + 2*X0*sqrt(1 + X0^2) )
      }
      if( length(param) == 1 & tick[1] == 0 )
      {
        X0 <- 1*rnorm(N*K)/2
        W[, g] <- theta[[1]][1]*( 1 + 2*X0^2 + 2*X0*sqrt(1 + X0^2) )
      }
      if( length(param) == 1 & tick[2] == 0 )
      {
        X0 <- theta[[g]][1]*rnorm(N*K)/2
        W[, g] <- 1*( 1 + 2*X0^2 + 2*X0*sqrt(1 + X0^2) )
      }
    }
    if(family == "chisq")  W[, g] <- rchisq(N*K, df = theta[[g]][1])
    if(family == "f")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- rf(N*K, df1 = theta[[g]][1], df2 = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- rf(N*K, df1 = theta[[g]][1], df2 = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- rf(N*K, df1 = 1            , df2 = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <- rf(N*K, df1 = theta[[g]][1], df2 = 1             )
    }
    if(family == "gamma")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- rgamma(N*K, shape = theta[[g]][1], rate = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- rgamma(N*K, shape = theta[[g]][1], rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- rgamma(N*K, shape = 1            , rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <- rgamma(N*K, shape = theta[[g]][1], rate = 1             )
    }
    if(family == "igamma")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- 1/rgamma(N*K, shape = theta[[g]][1], rate = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- 1/rgamma(N*K, shape = theta[[g]][1], rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- 1/rgamma(N*K, shape = 1            , rate = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <- 1/rgamma(N*K, shape = theta[[g]][1], rate = 1             )
    }
    if(family == "gigaussian")
    {
      if( length(param) == 3 & sum( tick ) == 3 ) W[, g] <- rgig( N*K, lambda = theta[[g]][1], chi = theta[[g]][2], psi = theta[[g]][3] )
      if( length(param) == 2 & sum( tick ) == 3 ) W[, g] <- rgig( N*K, lambda = theta[[g]][1], chi = theta[[g]][2], psi = theta[[g]][2] )
      if( length(param) == 2 &     tick[1] == 0 ) W[, g] <- rgig( N*K, lambda = 1, chi = theta[[g]][1], psi = theta[[g]][2] )
      if( length(param) == 2 &     tick[2] == 0 ) W[, g] <- rgig( N*K, lambda = theta[[g]][1], chi = 1, psi = theta[[g]][2] )
      if( length(param) == 2 &     tick[3] == 0 ) W[, g] <- rgig( N*K, lambda = theta[[g]][1], chi = theta[[g]][2], psi = 1 )
    }
    if(family == "igaussian")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- rigaussian( N*K, alpha = theta[[g]][1], beta = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- rigaussian( N*K, alpha = theta[[g]][1], beta = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- rigaussian( N*K, alpha = 1            , beta = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <- rigaussian( N*K, alpha = theta[[g]][1], beta = 1             )
    }
    if(family == "lidley")
    {
      X0 <- rbinom( N*K, 1, theta[[g]][1]/(1 + theta[[g]][1]) )
      W[, g] <- X0*rgamma( N*K, shape = 1, rate = theta[[g]][1] ) + (1 - X0)*rgamma( N*K, shape = 2, rate = theta[[g]][1] )
    }
    if(family == "loglog")
    {
      u <- runif( N*K )
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- theta[[g]][2]*( u/(1 - u) )^(1/theta[[g]][1] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- theta[[g]][1]*( u/(1 - u) )^(1/theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- theta[[g]][2]*( u/(1 - u) )^(1/1             )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <-             1*( u/(1 - u) )^(1/theta[[g]][1] )
    }
    if(family == "lognorm")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- rlnorm( N*K, meanlog = theta[[g]][1], sdlog = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- rlnorm( N*K, meanlog = theta[[g]][1], sdlog = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- rlnorm( N*K, meanlog = 1            , sdlog = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <- rlnorm( N*K, meanlog = theta[[g]][1], sdlog = 1             )
    }
    if(family == "lomax")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- 1/theta[[g]][2]*( (1 - runif( N*K ) )^(-1/theta[[g]][1]) - 1 )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- 1/theta[[g]][1]*( (1 - runif( N*K ) )^(-1/theta[[g]][1]) - 1 )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- 1/theta[[g]][1]*( (1 - runif( N*K ) )^(-1/1) - 1 )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <-             1/1*( (1 - runif( N*K ) )^(-1/theta[[g]][1]) - 1 )
    }
    if(family == "pstable")
    {
      w0     <- -log( runif(N*K) )
      theta0 <-  pi*( runif(N*K) - 1/2 )
      sigma0 <-  (cos(pi*theta[[g]][1]/4))^(2/theta[[g]][1])
      r0     <-  sin(theta[[g]][1]/2*(pi/2 + theta0))/( cos(theta[[g]][1]*pi/4)*cos(theta0) )^(2/theta[[g]][1])*
        ( cos(theta[[g]][1]*pi/4 + (theta[[g]][1]/2 - 1)*theta0)/w0 )^( (2 - theta[[g]][1])/theta[[g]][1] )
      W[, g] <-  ( r0 - tan(pi*theta[[g]][1]/4) )*sigma0 + sigma0*tan(pi*theta[[g]][1]/4)
    }
    if(family == "ptstable")  W[, g] <- 1/( 2^(1 - 2/theta[[g]][1])*rptstable(N*K, theta[[g]][1], Dim/2) )
    if(family == "rayleigh")  W[, g] <- rweibull(N*K, shape = 2, scale = theta[[g]][1] )
    if(family == "weibull")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- rweibull(N*K, shape = theta[[g]][1], scale = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- rweibull(N*K, shape = theta[[g]][1], scale = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- rweibull(N*K, shape = 1            , scale = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <- rweibull(N*K, shape = theta[[g]][1], scale = 1             )
    }
    if(family == "constant")
    {
      W[, g]    <- rep(1, N*4)
      T0[, , g] <- qnorm( ( runif(Q*N*4) + 1 )/2 )
      T1[, , g] <- sapply(1:(N*4), function(r) T0[r, , g]*sqrt( W[r, g] ) )
    }else{
      T0[, , g] <- qnorm( ( runif(Q*N*K) + 1 )/2 )
      T1[, , g] <- sapply(1:(N*K), function(r) T0[r, , g]*sqrt( W[r, g] ) )
    }
    ############ computing  f_y_g ########
    Mu     <- as.vector( mu[[g]] )
    Sigma  <- as.matrix( sigma[[g]] )
    Lambda <- as.matrix( lambda[[g]] )
    C0[g] <- 1/( (2*pi)^( Dim/2 )*sqrt( abs( det( Sigma ) ) ) )
    for(i in 1:n)
    {
      ############ computing f(y)                     ########
      #k <- 1
      index <- 1:N #sample( ( N*(k - 1) + 1 ):( k*N ), M)
      X6 <- T1[index, , g]
      X7 <- W[index, g]
      X1 <- t( sapply(1:M, function(r) Y[i, ] - Mu - Lambda%*%X6[r, ] ) )
      P1 <- exp( log( C0[g] ) + ( -Dim/2     )*log( X7 ) -0.5*mahalanobis( X1, rep(0, Dim ), Sigma )/X7 )
      f_y[i, g] <- mean( P1 )
      ############ end of computing f(y)              ########
      ############ computing E(W^(-1) | Y)            ########
     # k <- 2
     # index <- sample( ( N*(k - 1) + 1 ):( k*N ), M)
      X6 <- T1[index, , g]
      X7 <- W[index, g]
      X2 <- t( sapply(1:M, function(r) Y[i, ] - Mu - Lambda%*%X6[r, ] ) )
      P2 <- exp( log( C0[g] ) + ( -Dim/2 - 1 )*log( X7 ) -0.5*mahalanobis( X2, rep(0, Dim ), Sigma )/X7 )
      E_W1[i, g] <- mean( P2 )/f_y[i, g]
      ############ end of computing E(ZW^(-1) | Y)    ########
      ############ computing E(W^(-1)T | Y)           ########
     # k <- 3
     # index <- sample( ( N*(k - 1) + 1 ):( k*N ), M)
      X6 <- T1[index, , g]
      X7 <- W[index, g]
      X3 <- t( sapply(1:M, function(r) Y[i, ] - Mu - Lambda%*%X6[r, ] ) )
      P3 <- exp( log( C0[g] ) + ( -Dim/2 - 1 )*log( X7 ) -0.5*mahalanobis( X3, rep(0, Dim ), Sigma )/X7 )
      E_W1_T[i, , g] <- colMeans( X6*P3 )/f_y[i, g]
      #  E_W1_T[i, , g] <- colMeans( t(sapply(index, function(r)T1[r, , g]*P3[r] )), na.rm = TRUE )/f_y[i, g]
      #  E_W1_T[i, , g] <- apply( sapply(index, function(r) T1[r, , g]*P3[r] ), 1, mean, na.rm = TRUE )/f_y[i, g]
      ############ end of computing E(ZW^(-1)T | Y)   ########
      ############ computing E(W^(-1)T^tT | Y)        ########
    #  k <- 4
    #  index <- sample( ( N*(k - 1) + 1 ):( k*N ), M)
      X6 <- T1[index, , g]
      X7 <- W[index, g]
      X4 <- t( sapply(1:M, function(r) Y[i, ] - Mu - Lambda%*%X6[r, ] ) )
      P4 <- exp( log( C0[g] ) + ( -Dim/2 - 1 )*log( X7 ) -0.5*mahalanobis( X4, rep(0, Dim ), Sigma )/X7 )
      E_W1_TT[i, , g] <- rowMeans( sapply(1:M, function(r) X6[r, ]%o%X6[r, ]*P4[r] ) )/f_y[i, g]
      ############ end of computing E(W^(-1)T^tT | Y) ########
    }
  }
  ############ computing E(Z | Y)        ########
  tau_ig <- sweep( f_y, 2, weight, "*" )
  tau_ig <- tau_ig/rowSums( tau_ig )
  ############ end of computing E(Z | Y) ########
  for( g in 1:G )
  {
    Mu     <- as.vector( mu[[g]] )
    Sigma  <- as.matrix( sigma[[g]] )
    Lambda <- as.matrix( lambda[[g]] )
    E_W1[, g] <- E_W1[, g]*tau_ig[, g]
    E_W1_T[, , g] <- E_W1_T[, , g]*tau_ig[, g]
    E_W1_TT[, , g] <- E_W1_TT[, , g]*tau_ig[, g]
    ############ computing derivative of f_W(w) w.r.t theta ########
    if( any( family == c( "pstable", "ptstable" ) ) )
    {
    #  k <- 5
    #  index <- sample( ( N*(k - 1) + 1 ):( k*N ), M)
      X6 <- T1[index, , g]
      X7 <- W[index, g]
      if( family == "pstable")
    {
      forward_diff  <- dpstable( W[index, g], theta[[g]] + h )
      backward_diff <- dpstable( W[index, g], theta[[g]] - h )
      derive_theta  <- ( forward_diff - backward_diff )/( 2*h*dpstable( W[index, g], theta[[g]] ) )
    }else{
      forward_diff  <- dptstable( W[index, g], theta[[g]] + h )
      backward_diff <- dptstable( W[index, g], theta[[g]] - h )
      derive_theta  <- ( forward_diff - backward_diff )/( 2*h*dptstable( W[index, g], theta[[g]] ) )
    }
      for(i in 1:n)
      {
        X5 <- t( sapply(1:M, function(r) Y[i, ] - Mu - Lambda%*%X6[r, ] ) )
        P5 <- exp( log( C0[g] ) -( Dim/2 )*log( X7 ) - 0.5*mahalanobis( X5, rep(0, Q ), Sigma )/X7 )
        E_deriv_theta[i, 1, g] <- mean( derive_theta*P5, na.rm = TRUE)/f_y[i, g]
      }
    }
    if( all( family != c( "pstable", "ptstable", "constant" ) ) )
    {
      first_deriv <- tryCatch( sapply( 1:Dim_theta, function(r) D( bquote( .( PDF ) ), param[r] ) ), error = function(e)( "fail" )  )
      if( any( first_deriv == "fail" ) )
      {
        pdf1 <- function(w, param){ }
        body( pdf1 ) <- bquote( .( PDF ) )
        theta2[[g]] <- theta[[g]]
        theta1[[g]] <- theta[[g]]
        for(j in 1:Dim_theta)
        {
         # k <- 4 + j
         # index <- sample( ( N*(k - 1) + 1 ):( k*N ), M)
          X6 <- T1[index, , g]
          X7 <- W[index, g]
          theta2[[g]][j] <- theta[[g]][j] + h
          theta1[[g]][j] <- theta[[g]][j] - h
          for(r in 1:Dim_theta) assign( param[r], theta2[[g]][r] )
          forward_diff <- pdf1( X7, theta2[[g]] )
           for(r in 1:Dim_theta) assign( param[r], theta1[[g]][r] )
            backward_diff <- pdf1( X7, theta1[[g]] )
            for(r in 1:Dim_theta) assign( param[r], theta[[g]][r] )
            derive_theta <- ( forward_diff - backward_diff )/( 2*h*pdf1( X7, theta[[g]] ) )
          for(i in 1:n)
          {
            X5 <- t( sapply(1:M, function(r) Y[i, ] - Mu - Lambda%*%X6[r, ] ) )
            P5 <- exp( log( C0[g] ) -( Dim/2 )*log( X7 ) - 0.5*mahalanobis( X5, rep(0, Q ), Sigma )/X7 )
            E_deriv_theta[i, j, g] <- mean( derive_theta*P5, na.rm = TRUE)/f_y[i, g]
          }
          theta1[[g]] <- theta[[g]]
          theta2[[g]] <- theta[[g]]
        }
      }else{
        pdf2 <- function(w, y, i){ }
        for(r in 1:Dim_theta) assign( param[r], theta[[g]][r] )
        for(j in 1:Dim_theta)
        {
          body(pdf2) <- bquote(	.( first_deriv[[j]] )/.( PDF )*.( quote( exp( log( C0[g] ) -( Dim/2 )*log( w ) -
                                0.5*mahalanobis( y, rep(0, Q ), Sigma )/w ) ) )	)
          for(i in 1:n)
          {
          #  k <- 4 + j
          #  index <- sample( ( N*(k - 1) + 1 ):( k*N ), M)
            X6 <- T1[index, , g]
            X7 <- W[index, g]
            X5 <- t( sapply(1:M, function(r) Y[i, ] - Mu - Lambda%*%X6[r, ] ) )
            E_deriv_theta[i, j, g] <- mean( pdf2( X7, X5, i), na.rm = TRUE )/f_y[i, g]
          }
        }
      }
    }
    E_deriv_theta[, , g] <- E_deriv_theta[, , g]*tau_ig[, g]
  }
  ########### end of computing derivative of f_W(w) w.r.t theta ########
  for(g in 1:G)
  {
    for(i in 1:n)
    {
      if( is.nan( tau_ig[i, g] ) ) tau_ig[i, g]  <- 0
      if( is.nan( E_W1[i, g] ) | is.infinite( E_W1[i, g] ) )
      {
        E_W1[i, g] <- 0
      }
      for(j in 1:Dim_theta)
      {
        if( is.nan( E_deriv_theta[i, j, g] ) ) E_deriv_theta[i, j, g] <- 0
        if( is.infinite( E_deriv_theta[i, j, g] ) ) E_deriv_theta[i, j, g] <- 0
      }
      for(k in 1:Q)
      {
        if( is.nan( E_W1_T[i, k, g] ) | is.infinite( E_W1_T[i, k, g] ) )
        {
          E_W1_T[i, k, g] <- 0
        }
      }
      A <- matrix( E_W1_TT[i, , g], nrow = Q, ncol = Q )
      for(r in 1:Q)
      {
        for(s in 1:Q)
        {
          if( is.nan( A[r, s] ) | is.infinite( A[r, s] ) )
          {
            A[r, s] <- 0
          }
        }
      }
      E_W1_TT[i, , g] <- ( A + t(A) )/2
    }
  }
  return( list(tau_ig = tau_ig, E_W1 = E_W1, E_deriv_theta = E_deriv_theta, E_W1_TT = E_W1_TT,
               E_W1_T = E_W1_T) )
}
