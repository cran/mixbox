estep1 <- function(Y, G, weight, mu, sigma, lambda, family, skewness, param, theta, tick, h, N, PDF)
{
  Dim <- length( mu[[1]] )
  P1 <- C0 <- rep( 0, G )
  n <- length( Y[, 1] )
  M <- N #floor( 0.9*N )
  delta <- rep( 0, G )
  dy <- tau_ig <- m <- tau_ig <- E_W1 <- f_y <- matrix( NA, nrow = n, ncol = G )
  Sigmainv <- list( )
  theta01  <- theta02 <- vector( "list", G )
  if( family == "constant" )
  {
    Dim_theta <- 1
    K <- 4 + Dim_theta
    W  <- matrix( 1, nrow = N*4, ncol = G )
  }else{
    Dim_theta <- length( theta[[1]] )
    K <- 4 + Dim_theta
    X0 <- rep(NA, N*K)
    W  <- matrix( NA, nrow = N*K, ncol = G )
  }
  E_deriv_theta <- array( 0, c( n, Dim_theta, G ) )
  E_W1_Tg  <- E_1j_g <- array( 0, c( n, 2, G ) )
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
    if(family == "chisq") W[, g] <- rchisq(N*K, df = theta[[g]][1])
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
    if(family == "ptstable") W[, g] <- 1/( 2^(1 - 2/theta[[g]][1])*rptstable(N*K, theta[[g]][1], Dim/2) )
    if(family == "rayleigh") W[, g] <- rweibull(N*K, shape = theta[[g]][1], scale = theta[[g]][2] )
    if(family == "weibull")
    {
      if( length(param) == 2 & sum( tick ) == 2 ) W[, g] <- rweibull(N*K, shape = theta[[g]][1], scale = theta[[g]][2] )
      if( length(param) == 1 & sum( tick ) == 2 ) W[, g] <- rweibull(N*K, shape = theta[[g]][1], scale = theta[[g]][1] )
      if( length(param) == 1 &     tick[1] == 0 ) W[, g] <- rweibull(N*K, shape = 1            , scale = theta[[g]][1] )
      if( length(param) == 1 &     tick[2] == 0 ) W[, g] <- rweibull(N*K, shape = theta[[g]][1], scale = 1             )
    }
    ########### computing  f_y_g ########
    Mu     <- as.vector( mu[[g]] )
    Sigma  <- as.matrix ( sigma[[g]] )
    Lambda <- as.vector( lambda[[g]] )
    Omega  <- Sigma + Lambda%*%t( Lambda )
    delta[g] <- 1 - mahalanobis( Lambda, rep(0, Dim), Omega )
    Sigmainv[[g]] <- solve( Sigma )
    C0[g]   <- 2/sqrt( (2*pi)^( Dim + 1 )* abs( det( Sigma ) ) )
    dy[, g] <- mahalanobis( Y, Mu, Omega )
    m[, g]  <- sapply(1:n, function(i) as.numeric(t( Lambda )%*%solve(Omega)%*%c(Y[i, ] - Mu)) )
    s <- 0
    #k <- 1
    index <- 1:N #sample( ( N*(k - 1) + 1 ):( k*N ), M)
    P1[g] <- C0[g]*( 2*delta[g] )^( (0 + 1)/2 ) * gamma( (0 + 1)/2 )/2
    P2 <- colMeans( sapply(1:n, function(i){ W[index, g]^( (0 + 1)/2 - (Dim + 1)/2 + s)*exp( -0.5*dy[i, g]/W[index, g] )*(1 +
                                             pgamma( m[i, g]^2/(2*delta[g]*W[index, g]), shape = (0 + 1)/2, rate = 1)*(-1)^(0)*( sign(m[i, g]) )^(0 + 1))} ) )
    f_y[, g] <- P1[g]*P2
    ############ end of computing  f(y)   ########
    ############ computing E(ZW^(-1) | Y) ########
    s <- -1
    #k <- 2
    #index <- sample( ( N*(k - 1) + 1 ):( k*N ), M)
    P2 <- colMeans( sapply(1:n, function(i){ W[index, g]^( (0 + 1)/2 - (Dim + 1)/2 + s)*exp( -0.5*dy[i, g]/W[index, g] )*(1 +
                                             pgamma( m[i, g]^2/(2*delta[g]*W[index, g]), shape = (0 + 1)/2, rate = 1)*(-1)^(0)*( sign(m[i, g]) )^(0 + 1))} ) )
    E_W1[, g] <- P1[g]*P2/f_y[, g]
    ############ end of computing  E(ZW^(-1) | Y) ########
    ############ computing E(ZW^(-1)T | Y)        ########
    for( j in 1:2 )
    {
      #k <- j + 2
      #index <- sample( ( N*(k - 1) + 1 ):( k*N ), M)
      P3 <- C0[g]*( 2*delta[g] )^( (j + 1)/2 )*gamma( (j + 1)/2 )/2
      P2 <- colMeans( sapply(1:n, function(i){ W[index, g]^( (j + 1)/2 - (Dim + 1)/2 + s)*exp( -0.5*dy[i, g]/W[index, g] )*(1 +
                                               pgamma( m[i, g]^2/(2*delta[g]*W[index, g]), shape = (j + 1)/2, rate = 1)*(-1)^(j)*( sign(m[i, g]) )^(j + 1)) } ) )
      E_1j_g[, j,  g] <- P3*P2/f_y[, g]
    }
    E_W1_Tg[, 1, g] <- (m[, g]*E_W1[, g] + E_1j_g[, 1, g] )
    E_W1_Tg[, 2, g] <- (m[, g]^2*E_W1[, g] + 2*m[, g]*E_1j_g[, 1, g] + E_1j_g[, 2, g] )
  }
  ############ end of computing E(ZW^(-1)T | Y) ########
  ############ computing E(Z | Y)               ########
  tau_ig <- sweep( f_y, 2, weight, "*" )
  tau_ig <- tau_ig/rowSums( tau_ig )
  ############ end of computing E(Z | Y) ########
  for( g in 1:G )
  {
    E_W1[,       g] <-    E_W1[,   g]*tau_ig[, g]
    E_W1_Tg[, ,  g] <- E_W1_Tg[, , g]*tau_ig[, g]
    ############ computing derivative of f_W(w) w.r.t theta ########
    if( all( family != c( "pstable", "ptstable", "constant" ) ) )
    {
      first_deriv <- tryCatch( sapply( 1:Dim_theta, function(r) D( bquote( .( PDF ) ), param[r] ) ), error = function(e)( "fail" )  )
      pdf1 <- function(w, param){ }
      body( pdf1 ) <- bquote( .( PDF ) )
      if( any(first_deriv == "fail") )
      {
        theta01[[g]] <- theta[[g]]
        theta02[[g]] <- theta[[g]]
        for(j in 1:Dim_theta)
        {
          #k <- 4 + j
          #index <- sample( ( N*(k - 1) + 1 ):( k*N ), M)
          theta01[[g]][j] <- theta[[g]][j] - h
          theta02[[g]][j] <- theta[[g]][j] + h
            for(r in 1:Dim_theta) assign( param[r], theta02[[g]][r] )
            forward_diff <- pdf1( W[index, g], theta02[[g]] )
            for(r in 1:Dim_theta) assign( param[r], theta01[[g]][r] )
            backward_diff <- pdf1( W[index, g], theta01[[g]] )
            for(r in 1:Dim_theta) assign( param[r], theta[[g]][r] )
            derive_theta <- ( forward_diff - backward_diff )/( 2*h*pdf1( W[index, g], theta[[g]] ) )
            for(i in 1:n)
            {
            E_deriv_theta[i, j, g] <- 2*P1[g]*mean( derive_theta*W[index, g]^(-Dim/2)*exp(-dy[i, g]/(2*W[index, g]))*
                                                    pnorm( m[i, g] /sqrt( W[index, g]*delta[g] ) ), na.rm = TRUE)/f_y[i, g]
			      }
          theta01[[g]] <- theta[[g]]
          theta02[[g]] <- theta[[g]]
        }
      }else{
        pdf2 <- function(w, i){ }
        for(r in 1:Dim_theta) assign( param[r], theta[[g]][r] )
        for(j in 1:Dim_theta)
        {
          #k <- 4 + j
          #index <- sample( ( N*(k - 1) + 1 ):( k*N ), M)
          body(pdf2) <- bquote( .( first_deriv[[j]] )/.( PDF ) *.( quote( w^(-Dim/2) ) ) *.( quote( exp( -dy[i, g]/(2*w) )) ) *
                                .( quote( pnorm( ( m[i, g] )/sqrt(w*delta[g]) ) ) ) )
          for(i in 1:n)
          {
            E_deriv_theta[i, j, g] <- 2*P1[g]*mean( pdf2( W[index, g], i), na.rm = TRUE )/f_y[i, g]
          }
        }
      }
    }
    if( family == "pstable")
    {
        #k <- 5
        #index <- sample( ( N*(k - 1) + 1 ):( k*N ), M)
        forward_diff  <- dpstable( W[index, g], theta[[g]] + h )
        backward_diff <- dpstable( W[index, g], theta[[g]] - h )
        derive_theta  <- ( forward_diff - backward_diff )/( 2*h*dpstable( W[index, g], theta[[g]] ) )
        for(i in 1:n)
        {
          E_deriv_theta[i, , g] <- 2*P1[g]*mean( derive_theta*W[index, g]^(-Dim/2)*exp( -dy[i, g]/(2*W[index, g]) )*
                                                 pnorm( m[i, g] /sqrt( W[index, g]*delta[g] ) ), na.rm = TRUE)/f_y[i, g]
        }
    }
    if( family == "ptstable")
    {
      #k <- 5
      #index <- sample( ( N*(k - 1) + 1 ):( k*N ), M)
          forward_diff  <- dptstable( W[index, g], theta[[g]] + h )
          backward_diff <- dptstable( W[index, g], theta[[g]] - h )
          derive_theta  <- ( forward_diff - backward_diff )/( 2*h*dptstable( W[index, g], theta[[g]] ) )
          for(i in 1:n)
          {
          E_deriv_theta[i, , g] <- 2*P1[g]*mean( derive_theta*W[index, g]^(-Dim/2)*exp(-dy[i, g]/(2*W[index, g]))*
                                                 pnorm( m[i, g] /sqrt( W[index, g]*delta[g] ) ), na.rm = TRUE)/f_y[i, g]
          }
    }
    E_deriv_theta[, , g] <- E_deriv_theta[, , g]*tau_ig[, g]
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
      if( is.nan(tau_ig[i, g]) || is.nan(E_W1[i, g]) || is.nan(E_W1_Tg[i, 1, g]) || is.nan(E_W1_Tg[i, 2, g]) )
      {
        tau_ig[i, g] <- 0
        E_W1[i, g] <- 0
        E_W1_Tg[i, 1, g]<- 0
        E_W1_Tg[i, 2, g]<- 0
      }
    }
  }
  return( list(tau_ig = tau_ig, E_W1 = E_W1, E_deriv_theta = E_deriv_theta, E_W1_Tg = E_W1_Tg) )
}

