configuration1 <- function(Y, G, weight, mu, sigma, lambda, family, skewness, param, theta, ofim1_solve, sigma_arrange1, level)
{
  Dim <- length( mu[[1]] )
  Dim_sigma  <- Dim*(Dim + 1)/2
  Dim_theta  <- length( theta[[1]] )
  sd_omega   <- rep(NA, G - 1 )
  Z_alpha2   <- qnorm(1 - level/2)
  out_weight <- matrix( NA, nrow = G - 1    , ncol = 4 )
  out_mu     <- matrix( NA, nrow = Dim      , ncol = 4 )
  out_sigma  <- matrix( NA, nrow = Dim_sigma, ncol = 4 )
  out_lambda <- matrix( NA, nrow = Dim      , ncol = 4 )
  out_theta  <- matrix( NA, nrow = Dim_theta, ncol = 4 )
  out_Mu     <- vector("list", G)
  out_Sigma  <- vector("list", G)
  out_Lambda <- vector("list", G)
  out_Theta  <- vector("list", G)
  squre_root_FI <- sqrt( abs( diag( abs( ofim1_solve ) ) ) )# sqrt( diag( solve( ofim1$Fisher ) ) )
  for( g in 1:(G - 1) )
  {
    out_weight[g, ] <- cbind( weight[g], squre_root_FI[g], weight[g] - squre_root_FI[g]*Z_alpha2, weight[g] + squre_root_FI[g]*Z_alpha2 )
  }
  for( g in 1:G )
  {
    for(d in 1:Dim)
    {
      range_mu     <- seq( G + (g - 1)*Dim, G + g*Dim - 1 )
      range_sigma  <- seq( G + G*Dim + (g - 1)*Dim_sigma, G + G*Dim + g*Dim_sigma - 1 )
      range_lambda <- seq( G + G*Dim + G*Dim_sigma + (g - 1)*Dim, G + G*Dim + G*Dim_sigma + g*Dim - 1 )
      range_theta  <- seq( G + G*Dim + G*Dim_sigma + G*Dim + (g - 1)*Dim_theta, G + G*Dim + G*Dim_sigma + G*Dim + g*Dim_theta - 1 )
      out_mu[d, ]  <- cbind( mu[[g]][d],  squre_root_FI[ range_mu ][d],
                             mu[[g]][d] - squre_root_FI[ range_mu ][d]*Z_alpha2,
                             mu[[g]][d] + squre_root_FI[ range_mu ][d]*Z_alpha2
      )
    }
    seq_Dim_sigma <- rep(NA, Dim_sigma)
    index <- 1
    for(i in 1:Dim)
    {
      for(j in i:Dim)
      {
        seq_Dim_sigma[index] <- paste(i, sep="", j)
        index = index + 1
      }
    }
    for(d in 1:Dim_sigma)
    {
      index <- sigma_arrange1[d, ]# ofim1$index_sigma[d, ]
      out_sigma[d, ] <- cbind( sigma[[g]][ index[1], index[2] ],  squre_root_FI[ range_sigma ][d],
                               sigma[[g]][ index[1], index[2] ] - squre_root_FI[ range_sigma ][d]*Z_alpha2,
                               sigma[[g]][ index[1], index[2] ] + squre_root_FI[ range_sigma ][d]*Z_alpha2
      )
    }
    if( family == "constant" & skewness == "TRUE"){
      for(d in 1:Dim)
      {
        range_lambda <- seq( G + G*Dim + G*Dim_sigma + (g - 1)*Dim, G + G*Dim + G*Dim_sigma + g*Dim - 1 )
        out_lambda[d, ] <- cbind( lambda[[g]][d],  squre_root_FI[ range_lambda ][d],
                                  lambda[[g]][d] - squre_root_FI[ range_lambda ][d]*Z_alpha2,
                                  lambda[[g]][d] + squre_root_FI[ range_lambda ][d]*Z_alpha2
        )
      }
    }
    if( family != "constant" & skewness == "FALSE" )
    {
      for(d in 1:Dim_theta)
      {
        range_theta  <- seq( G + G*Dim + G*Dim_sigma + 0*Dim + (g - 1)*Dim_theta, G + G*Dim +
                               G*Dim_sigma + 0*Dim + g*Dim_theta - 1 )
        out_theta[d, ] <- cbind( theta[[g]][d],  squre_root_FI[ range_theta ][d],
                                 theta[[g]][d] - squre_root_FI[ range_theta ][d]*Z_alpha2,
                                 theta[[g]][d] + squre_root_FI[ range_theta ][d]*Z_alpha2
        )
      }
    }
    if( family != "constant" & skewness == "TRUE" )
    {
      for(d in 1:Dim)
      {
        out_lambda[d, ] <- cbind( lambda[[g]][d],  squre_root_FI[ range_lambda ][d],
                                  lambda[[g]][d] - squre_root_FI[ range_lambda ][d]*Z_alpha2,
                                  lambda[[g]][d] + squre_root_FI[ range_lambda ][d]*Z_alpha2
        )
      }
      for(d in 1:Dim_theta)
      {
        out_theta[d, ] <- cbind( theta[[g]][d],  squre_root_FI[ range_theta ][d],
                                 theta[[g]][d] - squre_root_FI[ range_theta ][d]*Z_alpha2,
                                 theta[[g]][d] + squre_root_FI[ range_theta ][d]*Z_alpha2
        )
      }
    }
    seq_Dim_theta <- seq( 1:Dim_theta)
    colnames( out_weight ) <- c("MLE", "SE", "lower bound", "upper bound")
    rownames( out_weight ) <- paste( expression(weight), sep = "", seq(1, G - 1 ) )
    colnames( out_mu )     <- c(paste( expression(MLE), sep = "" ), "SE", "lower bound", "upper bound")
    rownames( out_mu )     <-   paste( expression(mu), sep = "", g, 1:Dim   )
    out_Mu[[g]] <- out_mu
    colnames( out_sigma )  <- c(paste( expression(MLE), sep = "" ), "SE", "lower bound", "upper bound")
    rownames( out_sigma )  <-   paste( expression(sigma), sep = "", seq_Dim_sigma )
    out_Sigma[[g]] <- out_sigma
    if( family == "constant" & skewness == "TRUE" )
    {
      colnames( out_lambda ) <- c(paste( expression(MLE), sep = "" ), "SE", "lower bound", "upper bound")
      rownames( out_lambda ) <-   paste( expression(lambda), sep = "", g, 1:Dim )
      out_Lambda[[g]] <- out_lambda
    }
    if( family != "constant" & skewness == "FALSE" )
    {
      colnames( out_theta )  <- c(paste( expression(MLE), sep = "" ), "SE", "lower bound", "upper bound")
      rownames( out_theta )  <-  paste( expression(param), sep = "", g, seq_Dim_theta )
      out_Theta[[g]] <- out_theta
    }
    if( family != "constant" & skewness == "TRUE" )
    {
      colnames( out_lambda ) <- c(paste( expression(MLE)   , sep = "" ), "SE", "lower bound", "upper bound")
      rownames( out_lambda ) <-   paste( expression(lambda), sep = ""  , g, 1:Dim   )
      out_Lambda[[g]] <- out_lambda
      colnames( out_theta )  <- c(paste( expression(MLE)   , sep = "" ), "SE", "lower bound", "upper bound")
      rownames( out_theta )  <-   paste( (param)           , sep = "" )
      out_Theta[[g]] <- out_theta
    }
  }
  if( family == "constant" & skewness == "FALSE" ){out <- list(weight = out_weight, Mu = out_Mu, Sigma = out_Sigma )}
  if( family == "constant" & skewness == "TRUE"  ){out <- list(weight = out_weight, Mu = out_Mu, Sigma = out_Sigma, Lambda = out_Lambda )}
  if( family != "constant" & skewness == "FALSE" ){out <- list(weight = out_weight, Mu = out_Mu, Sigma = out_Sigma, Theta = out_Theta )}
  if( family != "constant" & skewness == "TRUE"  ){out <- list(weight = out_weight, Mu = out_Mu, Sigma = out_Sigma, Lambda = out_Lambda, Theta = out_Theta )}
  return(out)
}

