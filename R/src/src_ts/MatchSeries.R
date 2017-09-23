#-------------------------------------------------------------------------------
MatchSeries <- function(t1, t2, smth.n = 5){
  
  t1.tr <- range(t1[, 1])
  t1.stp <- GetSamplingIntervalMode(t1[, 1], dg = 5)
  
  t2.tr <- range(t2[, 1])
  t2.stp <- GetSamplingIntervalMode(t2[, 1], dg = 5)
  
  
  # Establish temporal window where analysis make sense
  # Cut based on maximum minimum
  min.t <- max(c(t1.tr[1], t2.tr[1]))
  max.t <- min(c(t1.tr[2], t2.tr[2]))
  t.stp <- max(t1.stp, t2.stp)
  
  # Test if time steps are all the smae
  if((t.stp == t1.stp) & (t.stp == t2.stp)){
    if(!any(is.na(match(t2[, 1], t1[, 1])))){
      # Test which series is shorter
      if(length(t1[, 1]) >= length(t2[, 1])){
        out.indxs <- match(t2[, 1], t1[, 1])
        t1.out <- t1[out.indxs, ]
        t2.out <- t2
      } else{
        out.indxs <- match(t1[, 1], t2[, 1])
        t2.out <- t2[out.indxs, ]
        t1.out <- t1
      }
    } else{
      # Time to interpolate
      t.interpol <- seq(from = min.t, to = max.t, by = t.stp)
      # Interpolate
      t1.out <- InterpolateWithLoess(dat = t1, smth.n = smth.n,
                                     interpol.times = t.interpol)
      t2.out <- InterpolateWithLoess(dat = t2, smth.n = smth.n,
                                     interpol.times = t.interpol)
    }
    
  } else{
    # Time to interpolate
    t.interpol <- seq(from = min.t, to = max.t, by = t.stp)
    # Interpolate
    t1.out <- InterpolateWithLoess(dat = t1, smth.n = smth.n,
                                   interpol.times = t.interpol)
    t2.out <- InterpolateWithLoess(dat = t2, smth.n = smth.n,
                                   interpol.times = t.interpol)
  }
  
  #as.data.frame(approx(t1, xout = t.interpol, method = "linear"))
  #plot(t1,type="o")
  #points(t1.out,col="red",pch=19,cex=0.5)
  #as.data.frame(approx(t2, xout = t.interpol, method = "linear"))
  #plot(t2,type="o",ylim=c(-0.25,0.25))
  #points(t2.out,col="red",pch=19,cex=0.5,type="o")
  
  return(list(t1.out, t2.out))
}
#-------------------------------------------------------------------------------