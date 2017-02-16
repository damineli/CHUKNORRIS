#-----------------------------------------------------------------------------
GetSmoothVec <- function(sampling.freq, series.length, low.per = 16, 
                         high.per = 256){
  #in seconds!
  per.bands <- GetFreqOrPerBands(sampling.freq, 
                                 band.n = 11, 
                                 per = TRUE)
  
  per_band_max <- per.bands$per.band.max
  per_band_min <- per.bands$per.band.min
  
  vec.length <- sum(per_band_max <= series.length)
  smth.vec <- rep(TRUE, vec.length)
  
  smth.vec[(per_band_max[1:vec.length] <= low.per)] <- FALSE
  smth.vec[(per_band_min[1:vec.length] >= high.per)] <- FALSE
  
  return(smth.vec)
}
#-----------------------------------------------------------------------------