#-------------------------------------------------------------------------------
PreFilterTimeSeries <- function(dat, ts.par, plot.par){
  # Temporarily transfer time scale to seconds
  #if(ts.par$time.unit == "minutes"){dat[, 1] <- dat[, 1] * 60}
  # Intepolate if time steps are not homogeneous
  if(TestSamplingIntervalHomogeneity(dat[, 1], dg = ts.par$time.dgts)){
    time.round <- round(dat[, 1], digits = ts.par$time.dgts)
    dat.prefilt <- cbind(time.round, dat[, 2])
    colnames(dat.prefilt) <- colnames(dat)
    # Establish time step
    time.step <- dat.prefilt[2, 1] - dat.prefilt[1, 1]
  } else{
    # Establish time step
    time.step <- GetSamplingIntervalMode(dat[, 1], dg = ts.par$time.dgts)
    # Determine points where to interpolate
    interpol.times <- seq(from = min(dat[,1],na.rm = T), 
                          to = max(dat[,1],na.rm = T),
                          by = time.step)
    # Interpolate / lightly smooth ###PLOTTING POINT!
    dat.prefilt <- InterpolateWithLoess(dat,
                                        smth.n = ts.par$interpol.n,
                                        interpol.times = interpol.times)
  }
  if(ts.par$outlier.removal != "none"){
    # Remove outliers ###PLOTTING POINT!
    dat.clean <- RemoveOutliers(dat.prefilt, sampling.freq = 1 / time.step,
                                outlier.removal = ts.par$outlier.removal, 
                                integrate.ts = ts.par$outl.cum.ts, 
                                out.thrsh = ts.par$out.thrsh,
                                smth.n = ts.par$interpol.n)
    
  } else{dat.clean <- dat.prefilt}
  if(ts.par$coarse.smth){
    def.spn <- ts.par$coarse.n / dim(dat.clean)[1]
    y.new <- IntegrateAndSmoothWithLoess(dat.clean[, 2],
                                         spn =  def.spn,
                                         only.positive = ts.par$coarse.positive,
                                         integrate.ts = ts.par$coarse.integrate)
    dat.clean <- cbind(dat.clean[, 1], y.new)
  }
  
  colnames(dat.clean) <- colnames(dat)
  
  if(plot.par$do.plot){
    if(ts.par$time.unit == "min"){plt.fac <- 1 / 60} else{plt.fac <- 1}
    plot(dat[ ,1] * plt.fac, dat[ ,2], type="o", col = "grey", lwd = 1.5, 
         xlab = paste("Time (", ts.par$time.unit, ")", sep = ""), 
         ylab = plot.par$ylab, main = "Prefiltering results")
    points(dat.prefilt[ ,1] * plt.fac, dat.prefilt[ ,2], 
           col = adjustcolor("red", alpha.f = 0.75), pch = 19, cex = 0.5,
           type = "o")
    points(dat.clean[ ,1] * plt.fac, dat.clean[ ,2], 
           col = adjustcolor("green", alpha.f = 0.75), pch = 19, cex = 0.5, 
           type = "o")
    legend("bottomleft", legend = c("Raw", "Interpolated", "Smoothed"), 
           col = c("grey", "red","green"), bty = "n", lwd = c(1.5, 1, 1), 
           pch = c(1, 19, 19), pt.cex = c(1, 0.5, 0.5))
    }
  
  return(dat.clean)
}
#-------------------------------------------------------------------------------