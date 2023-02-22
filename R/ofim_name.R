ofim_name <- function(G, weight, Dim, lambda, model = "restricted", family = "constant", skewness = "FALSE", param )
{
  Dim_sigma  <- Dim*(Dim + 1)/2
  Dim_lambda <- Dim
  sd_omega   <- rep(NA, G - 1 )
  if( family != "constant" )
  {
    Dim_theta <- length( param )
  }else{
    Dim_theta <- 0
  }
  if(model == "canonical" & skewness == "TRUE" )
  {
    Q <- length( lambda[[1]][1,] )
    Dim_lambda <- Dim*Q
  }else{
    if(model != "canonical" & skewness == "TRUE" )
    {
      Dim_lambda <- Dim
    }else{
      Dim_lambda <- 0
    }
  }
  S <- rep( 0, G - 1 + G*( Dim + Dim_lambda + Dim_sigma + Dim_theta ) )
  S_name <- rep(NA, length(S) )
  for(g in 1:G)
  {
    seq_Dim_sigma <- rep(NA, Dim_sigma)
    index <- 1
    for(i in 1:Dim)
    {
      for(j in i:Dim)
      {
        seq_Dim_sigma[index] <- paste(g,".",i, sep="", j)
        index <- index + 1
      }
    }
    if(model == "canonical" & skewness == "TRUE")
    {
      seq_Dim_lambda <- rep(NA, Dim_lambda)
      index <- 1
      for(i in 1:Dim)
      {
        for(j in 1:Q)
        {
          seq_Dim_lambda[index] <- paste(g,".",i, sep="", j)
          index <- index + 1
        }
      }
    }
    if(model != "canonical" & skewness == "TRUE") seq_Dim_lambda <- 1:Dim
    range_mu     <- seq( G + (g - 1)*Dim, G + g*Dim - 1 )
    range_sigma  <- seq( G + G*Dim + (g - 1)*Dim_sigma, G + G*Dim + g*Dim_sigma - 1 )
    range_lambda <- seq( G + G*Dim + G*Dim_sigma + (g - 1)*Dim_lambda, G + G*Dim +
                           G*Dim_sigma + g*Dim_lambda - 1 )
    range_theta  <- seq( G + G*Dim + G*Dim_sigma + G*Dim_lambda + (g - 1)*Dim_theta,
                         G + G*Dim + G*Dim_sigma + G*Dim_lambda + g*Dim_theta - 1  )
    if(g < G ) S_name[ g ] <- paste( expression(weight), sep = "", g )
    S_name[ range_mu     ] <- paste( expression(mu    ), sep = "", g, 1:Dim       )
    S_name[ range_sigma  ] <- paste( expression(sigma ), sep = "", seq_Dim_sigma  )
    if(model == "canonical" & skewness == "TRUE") S_name[ range_lambda ] <- paste(
      expression(lambda), sep = "" , seq_Dim_lambda )
    if(model != "canonical" & skewness == "TRUE") S_name[ range_lambda ] <- paste(
      expression(lambda), sep = "" , g, 1:Dim       )
    if( family != "constant" ) S_name[ range_theta ] <- paste((param), sep = "", g)
  }
  return(S_name)
}
