#-------------------------------------------------------------------------------
RobustGMM <- function(obs, mu = NULL, verb = FALSE){
  out.gmm <- NULL
  if(length(mu) == 1){mu <- NULL}
  try(capture.output(out.gmm <- normalmixEM(x = obs, mu = mu, 
                             maxit = 2000, maxrestarts = 250, verb = verb))) 
if(!is.null(out.gmm)){
  try.n <- 1
  max.try <- 5
  while(is.null(out.gmm) & (try.n < max.try)){
    if(verb){message(paste("Attempt", try.n))}
    try(capture.output(out.gmm <- normalmixEM(x = obs, mu = mu, 
                               maxit = 2000, maxrestarts = 250, verb = verb)))
        #RunGMM(obs, mu = mu, maxit = 2000, maxrestarts = 250,
    try.n <- try.n + 1
  } 
}
if(verb){message(paste("I've tried",try.n,"times!"))}
#if(try.n >= max.try){out.gmm <- NULL}
return(out.gmm)  
}
#-------------------------------------------------------------------------------