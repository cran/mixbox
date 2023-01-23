welcome<- function(){
  msg <- c(paste0(
    "Welcome to
 _____             _____     ____    ____               ____     ____                                          ____               ____
|     |           |     |   |    |  |    |             |    |   |    |                                        |    |             |    |
|      |         |      |   |____|   |    |           |    |    |    |                                         |    |           |    |
|     _ |       | _     |    ____     |    |         |    |     |    |                                          |    |         |    |
|    | | |     | | |    |   |    |     |    |       |    |      |    |                                           |    |       |    |
|    |  | |___| |  |    |   |    |      |    |     |    |       |    |                  __________________        |    |     |    |
|    |   |     |   |    |   |    |       |    |___|    |        |    |_____________    |                  |        |    |___|    |
|    |    |___|    |    |   |    |        |    ___    |         |                  |   |     ________     |         |    ___    |
|    |             |    |   |    |       |    |   |    |        |     ________     |   |    |        |    |        |    |   |    |
|    |             |    |   |    |      |    |     |    |       |    |        |    |   |    |        |    |       |    |     |    |
|    |             |    |   |    |     |    |       |    |      |    |        |    |   |    |        |    |      |    |       |    |
|    |             |    |   |    |    |    |         |    |     |    |________|    |   |    |________|    |     |    |         |    |
|    |             |    |   |    |   |    |           |    |    |                  |   |                  |    |    |           |    |
|____|             |____|   |____|  |____|             |____|   |__________________|   |__________________|   |____|             |____|version ",

    packageVersion("mixbox")),"\nType 'citation(\"mixbox\")' for citing this R package in publications.")
  return(msg)
}
.onAttach <- function(libname, pkgname) {
  mess <- welcome()
  if(!interactive())
    mess[1] <- paste("Package 'mixbox' version", packageVersion("mixbox"))
  packageStartupMessage(mess)
  invisible()
}
