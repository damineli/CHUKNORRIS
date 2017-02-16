################################################################################
#-------------------------------------------------------------------------------
IntegrateAndSmoothWithLoess <- function(y, spn = 0.03, degree = 2, 
                                        only.positive = FALSE, 
                                        integrate.ts = TRUE){
  if(integrate.ts){
    if(only.positive){
      y.min <- min(y, na.rm = TRUE)
      y <- y - y.min
    }
    y.in <- cumsum(y)
  } else{y.in <- y}
  
  y.loess <- SmoothWithLoess(y.in, spn, degree)
  
  if(integrate.ts){
    new.y <- c(y.loess[1], diff(y.loess))
    if(only.positive){new.y <- new.y + y.min}
    } else{new.y <- y.loess}
  
  return(new.y)
}
#-------------------------------------------------------------------------------
SmoothWithLoess <- function(y, spn = 0.03, degree = 2){
  x <- 1:length(y)
  smth <- loess(y ~ x, span = spn, control = loess.control(surface = "direct"), 
                family = "gaussian", degree = degree) # span controls smoothing
  #to increase accuracy, pass further arguments to loess.control,
  #statistics = "exact", trace.hat = "exact", iterations = 10
  return(smth$fitted)
}
#-------------------------------------------------------------------------------
ApplyLoess <- function(dat, spn = 0.03, degree = 2){
  return(loess(dat[, 2] ~ dat[, 1],span = spn,
               control = loess.control(surface = "direct"),
               family = "gaussian", degree = degree)) 
  #to increase accuracy, pass further arguments to loess.control,
  #statistics = "exact", trace.hat = "exact", iterations = 10
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Yes, I know. These are a bit too redundant but hey, it means it can be easily improved
#-------------------------------------------------------------------------------
InterpolateWithLoess <- function(dat, smth.n = 5, interpol.times = NULL){
  # From number of points to time series fraction
  def.spn <- smth.n / dim(dat)[1]
  # Calculate time step, if not supplied
  if(is.null(interpol.times)){
    time.step <- GetSamplingIntervalMode(dat[, 1])
    # Generate sequence of time points where to approximate the time series
    interpol.times <- seq(from = min(dat[,1],na.rm = T), 
                          to = max(dat[,1],na.rm = T),
                          by = time.step)
  }
  
  # Perform loess fit
  loess.mdl <- ApplyLoess(dat, spn = def.spn)
  # Approximate to interpolation times
  mdl.interp <- predict(loess.mdl, interpol.times, se = F)
  
  new.dat <- cbind(interpol.times, mdl.interp)
  colnames(new.dat) <- colnames(dat)
  
  return(new.dat)
}

################################################################################
