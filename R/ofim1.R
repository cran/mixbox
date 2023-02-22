ofim1 <- function(Y, G, weight, mu, sigma, lambda, family = "constant", skewness = "FALSE", param = NULL, theta = NULL, tick = NULL, h = 0.001, N = 3000, level = 0.05, PDF = NULL)
{
  n   <- length( Y[, 1] )
  Dim <- length( mu[[1]] )
  Dim_sigma <- Dim*(Dim + 1)/2
  if( family == "constant" )
  {
    Dim_theta <- 1
  }else{
    Dim_theta <- length( theta[[1]] )
  }
  if( family != "constant" & skewness ==  "TRUE" ) S <- rep(0, G - 1 + G*( Dim +   Dim + Dim_sigma +   Dim_theta ) )
  if( family == "constant" & skewness == "FALSE" ) S <- rep(0, G - 1 + G*( Dim +     0 + Dim_sigma + 0*Dim_theta ) )
  if( family != "constant" & skewness == "FALSE" ) S <- rep(0, G - 1 + G*( Dim +     0 + Dim_sigma +   Dim_theta ) )
  if( family == "constant" & skewness ==  "TRUE" ) S <- rep(0, G - 1 + G*( Dim +   Dim + Dim_sigma + 0*Dim_theta ) )
  z_alpha2 <- qnorm(1 - level/2)
  Fisher   <- matrix( 0, nrow = length(S), ncol = length(S) )
  Sigmainverse  <- array( NA, c(Dim, Dim, G) )
  E_W1_Tg 	    <- array( 0, c(n, 2, G) )
  E_deriv_theta <- array( 0, c(n, Dim_theta, G) )
  tau_ig  <- E_W1 <- matrix( 0, nrow = n, ncol = G )
  estep   <- estep1(Y = Y, G = G, weight = weight, mu = mu, sigma = sigma, lambda = lambda, family = family, skewness = skewness,
                 param = param, theta = theta, tick = tick, h = h, N = N, PDF = PDF)
  tau_ig  <- estep$tau_ig
  E_W1    <- estep$E_W1
  E_W1_Tg <- estep$E_W1_Tg
  E_deriv_theta <- estep$E_deriv_theta
  for(g in 1:G) Sigmainverse[, , g] <- solve( sigma[[g]] )
  for(g in 1:G)
  {
    range_mu     <- seq( G + (g - 1)*Dim, G + g*Dim - 1 )
    range_sigma  <- seq( G + G*Dim + (g - 1)*Dim_sigma, G + G*Dim + g*Dim_sigma - 1 )
    range_lambda <- seq( G + G*Dim + G*Dim_sigma + (g - 1)*Dim, G + G*Dim + G*Dim_sigma + g*Dim - 1 )
    range_theta  <- seq( G + G*Dim + G*Dim_sigma + G*Dim + (g - 1)*Dim_theta, G + G*Dim + G*Dim_sigma + G*Dim + g*Dim_theta - 1 )
    Sigmainv <- Sigmainverse[, , g]
    Mu <- as.vector( mu[[g]] )
    Lambda <- as.vector( lambda[[g]] )
    for(i in 1:n)
    {
      S1 <- E_W1[i,    g]*( c(Y[i, ] - Mu)%o%c(Y[i, ] - Mu) ) + E_W1_Tg[i, 2, g]*( Lambda%o%Lambda ) -
            E_W1_Tg[i, 1, g]*(Lambda%o%c(Y[i, ] - Mu) )  - E_W1_Tg[i, 1, g]*(c(Y[i, ] - Mu )%o%Lambda )
      S6 <- Sigmainv%*%( E_W1[i,    g]*c(Y[i, ] - Mu) - E_W1_Tg[i, 1, g]*Lambda )
      S7 <- Sigmainv%*%( E_W1_Tg[i, 1, g]*c(Y[i, ] - Mu) - E_W1_Tg[i, 2, g]*Lambda )
      out_sigma <- arrange_sigma( 1/2*( -tau_ig[i, g]*Sigmainv + Sigmainv%*%matrix( S1, nrow = Dim, ncol = Dim )%*%Sigmainv ) )
      S[ range_sigma ] <- out_sigma$y # S_Sigma
      if( family != "constant" & skewness == "TRUE" )
      {
        S[ range_lambda ] <- as.vector( S7 ) # S_Lambda
        outE_deriv_theta <- estep$E_deriv_theta[i, , g]
        outE_deriv_theta[ is.na( outE_deriv_theta ) ] <- 0
        S[ range_theta ] <- outE_deriv_theta # S_Theta
      }
      if( family != "constant" & skewness == "FALSE" )
      {
        range_theta <- seq( G + G*Dim + G*Dim_sigma + (g - 1)*Dim_theta, G + G*Dim + G*Dim_sigma + g*Dim_theta - 1 )
        S[ range_theta ] <- estep$E_deriv_theta[i, , g]
      }
      if( family == "constant" & skewness == "TRUE"  ) S[ range_lambda ] <- as.vector( S7 )
      if( family == "constant" & skewness == "FALSE" ) S <- S[seq( 1 : ( G - 1 + G*Dim + G*Dim_sigma ) )]
      if(g < G ) S[ g ] <- tau_ig[i, g]/weight[[g]] - tau_ig[i, G]/weight[[G]] # S_weight
      S[ range_mu ] <- as.vector( S6 ) # S_Mu
      Fisher <- Fisher + S%o%S
    }
  }
  Name <- ofim_name(G, weight, Dim, lambda, model = "restricted", family, skewness, param)
  colnames( Fisher ) <- Name
  rownames( Fisher ) <- Name
  return( list ( OFI = Fisher, index_sigma = out_sigma$index ) )
}


