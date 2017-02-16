#-------------------------------------------------------------------------------
GetFreqOrPerBands <- function(sampling.freq, band.n = 11, per = FALSE){
  
  freq.band.min <- sampling.freq / (2^(1:band.n + 1)) #2^(1:band.n))
  freq.band.max <- sampling.freq / (2^(1:(band.n))) #(2^(0:(band.n - 1)))
  approx.freq.band.max <- sampling.freq / (2^(1:band.n + 1)) ##2^(1:band.n))
  
  out.bands <- cbind.data.frame(freq.band.max, 
                                freq.band.min, 
                                approx.freq.band.max)
  
  if(per){
    out.bands <- 1 / out.bands
    colnames(out.bands) <- c("per.band.min", 
                             "per.band.max", 
                             "approx.per.band.min")
  }
  
  return(out.bands)
}
#-----------------------------------------------------------------------------