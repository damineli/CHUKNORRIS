#-------------------------------------------------------------------------------
DetectTipAndGrowthRate <- function(kymo, tip.find.par, image.par){
  time.vec <- 0:(dim(kymo)[1] - 1) * image.par$time.step
  
  # Get tip location series
  tip.loc <- FindTip(kymo = kymo,
                     kymo.span.n = tip.find.par$kymo.span.n,
                     kymo.loess.dg = tip.find.par$kymo.loess.dg,
                     tip.find.n = tip.find.par$tip.find.n,
                     fluo.chunk.n = tip.find.par$fluo.chunk.n,
                     fluo.thrsh.frac = tip.find.par$fluo.thrsh.frac,
                     tip.span.n = tip.find.par$tip.span.n,
                     tip.loess.dg = tip.find.par$tip.loess.dg,
                     use.smooth = tip.find.par$use.smooth,
                     fix.slope = tip.find.par$fix.slope,
                     rm.tip.out = tip.find.par$rm.tip.out,
                     rm.tip.out.spn = tip.find.par$rm.tip.out.spn,
                     rm.tip.out.dg = tip.find.par$rm.tip.out.dg)
  
  # Remove outliers before proceeding
  # if(tip.find.par$rmv.tip.out){
  	# y.old <- cumsum(tip.loc[, 2])
    # tso.mdl <- tsoutliers::tso(ts(y.old, frequency =  1/((time.vec[2]-time.vec[1])*60)), 
                               # types = c("LS", "TC")) # to try harder add:
                                # # ",maxit = 100, maxit.iloop = 10" to try harder
    # y.new <- tso.mdl$yadj
    # y.new <- c(y.new[1], diff(y.new))
    # dim(y.new)
    # tip.loc[, 2] <- y.new

  # }
  # Calculate growth rate time series for raw and smoothed tip location series
  growth.rate.raw <- CalculateGrowthRate(tip.loc = tip.loc$raw,
                                         pixel.size = image.par$pixel.size, 
                                         time.step = image.par$time.step)
  
  growth.rate.smth <- CalculateGrowthRate(tip.loc = tip.loc$smth,
                                          pixel.size = image.par$pixel.size,
                                          time.step = image.par$time.step)
  
  tip.dat <- cbind.data.frame("time" = time.vec, 
                              "tip.loc.raw" = tip.loc$raw, 
                              "tip.loc.smth" = tip.loc$smth,
                              "growth.raw" = c(NA, growth.rate.raw),
                              "growth.smth" = c(NA, growth.rate.smth))
  
  return(tip.dat)
}
#-------------------------------------------------------------------------------
FindTip <- function(kymo, kymo.span.n = 7, kymo.loess.dg = 2, tip.find.n = 9, 
                    fluo.chunk.n = 9, fluo.thrsh.frac = 0.4, tip.span.n = 7, 
                    tip.loess.dg = 2, use.smooth = TRUE, fix.slope = FALSE, 
                    rm.tip.out = FALSE, rm.tip.out.spn = 0.25, rm.tip.out.dg = 2){
  # Smooth image and return to the original matrix orientation
  if(!is.null(use.smooth)){
    if(use.smooth){
      imaj <- t(apply(kymo, 1, SmoothWithLoess, #INTEGRATE VERSION?
                      spn = kymo.span.n / (dim(kymo)[1]),
                      degree = kymo.loess.dg))
    } else{imaj <- kymo}
  } else{imaj <- kymo}
  
  imaj <- imaj - min(imaj)
  
  if(!is.null(fix.slope)){
    if(fix.slope){
      median_slope <- median(unlist(lapply(apply(imaj, 1, FindTipAt, tip.find.n, fluo.chunk.n, 
                                                 fluo.thrsh.frac, return.indxs = T), `[[`, 3)))
      
      tip.loc <- apply(imaj, 1, FindTipWithSlope, slope = median_slope,
                       n.pts = tip.find.n, min.chunk.size = fluo.chunk.n, 
                       fluo.thrsh.frac = fluo.thrsh.frac)
    } else{
      tip.loc <- apply(imaj, 1, FindTipAt, n.pts = tip.find.n, 
                       min.chunk.size = fluo.chunk.n, 
                       fluo.thrsh.frac = fluo.thrsh.frac)
    }
  } else{
    tip.loc <- apply(imaj, 1, FindTipAt, n.pts = tip.find.n, 
                     min.chunk.size = fluo.chunk.n, 
                     fluo.thrsh.frac = fluo.thrsh.frac)
  }

  if(rm.tip.out){
  	tip.loc <- RemoveTipOutliers(tip.loc, dg = rm.tip.out.dg, spn = rm.tip.out.spn, tol = 3.5)
  }
  tip.loc.smth <- SmoothWithLoess(y = tip.loc,
                                  spn = tip.span.n / length(tip.loc),
                                  degree = tip.loess.dg)
    
  
    # IntegrateAndSmoothWithLoess(y = tip.loc,
    #                                           spn = tip.span.n / length(tip.loc),
    #                                           degree = tip.loess.dg, 
    #                                           only.positive = FALSE, 
    #                                           integrate.ts = TRUE)
  return(cbind.data.frame("raw" = tip.loc, "smth" = tip.loc.smth))
}
#-------------------------------------------------------------------------------
FindTipAt <- function(t.slice, n.pts = 11, min.chunk.size = 10, 
                      fluo.thrsh.frac = 0.4, return.indxs = FALSE){
  
  fluo.diff <- diff(t.slice)
  
  inis <- which(diff(fluo.diff < 0) == 1) + 1
  ends <- which(diff(fluo.diff < 0) == -1) + 1 
  #plot(t.slice)
  #abline(v=inis,col="green")
  #abline(v=ends,col="red")
  
  if((length(inis) > 0) & (length(ends) > 0)){
    
    if(inis[length(inis)] > ends[length(ends)]){
      ends <- c(ends, length(fluo.diff))
    }
    if(inis[1] > ends[1]){
      inis <- c(1, inis)
    }
    
    chunk.sizes <- ends - inis
    
    chunk.max.fluo <- c()
    for(i in 1:length(inis)){
      chunk.max.fluo <- c(chunk.max.fluo, max(t.slice[inis[i]:ends[i]]))
    }
    
    big.drops <- which((chunk.sizes > min.chunk.size) & (chunk.max.fluo > (fluo.thrsh.frac * max(t.slice))))
    
    if(length(big.drops) > 0){
      big.drops.indx <- max(big.drops)
      tip.ini <- inis[big.drops.indx]
      tip.end <- ends[big.drops.indx]
      #abline(v=c(tip.ini,tip.end),lwd=4)
    } else{
      warning("No big chunks of fluorescence decrease!")
      tip.ini <- tip.end <- which.min(fluo.diff)
    }
  } else{
    warning("No fluorescence decreases on the time slice!")
    tip.ini <- tip.end <- which.min(fluo.diff)
  }
  
  fluo.indxs <- tip.ini:tip.end
  
  if(length(fluo.indxs) > 1){
    sharpest.slope <- fluo.indxs[which.min(diff(t.slice[fluo.indxs]))]
  } else{sharpest.slope <- fluo.indxs}
  #SLOW VERSION - FIND MIN INTERCEPT IN MULTIPLE FITS
  fit.vec <- (sharpest.slope - n.pts):(sharpest.slope + 1)
  
  slope.sharpness <- c()
  for(i in fit.vec){
    slope.sharpness <-c(slope.sharpness, 
                        sum(diff(t.slice[i:(i + n.pts)])))
  }
  
  fit.indxs <- fit.vec[which.min(slope.sharpness)]:(fit.vec[which.min(slope.sharpness)] + (n.pts - 1))
  mdl <- lm(fit.indxs ~  t.slice[fit.indxs])
  
  if(return.indxs){
    return(list("indxs" = fit.indxs, "intercept" = mdl$coefficients[1], "slope" = mdl$coefficients[2]))
  }else{return(mdl$coefficients[1])}
}
#-------------------------------------------------------------------------------
CalculateGrowthRate <- function(tip.loc, pixel.size, time.step){
  return((diff(tip.loc) * pixel.size) / time.step)
} 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
FindTipWithSlope <- function(t.slice, slope, n.pts = 11, min.chunk.size = 10, 
                      fluo.thrsh.frac = 0.4){
  
  fluo.diff <- diff(t.slice)
  
  inis <- which(diff(fluo.diff < 0) == 1) + 1
  ends <- which(diff(fluo.diff < 0) == -1) + 1 
  #plot(t.slice)
  #abline(v=inis,col="green")
  #abline(v=ends,col="red")
  
  if((length(inis) > 0) & (length(ends) > 0)){
    
    if(inis[length(inis)] > ends[length(ends)]){
      ends <- c(ends, length(fluo.diff))
    }
    if(inis[1] > ends[1]){
      inis <- c(1, inis)
    }
    
    chunk.sizes <- ends - inis
    
    chunk.max.fluo <- c()
    for(i in 1:length(inis)){
      chunk.max.fluo <- c(chunk.max.fluo, max(t.slice[inis[i]:ends[i]]))
    }
    
    big.drops <- which((chunk.sizes > min.chunk.size) & (chunk.max.fluo > (fluo.thrsh.frac * max(t.slice))))
    
    if(length(big.drops) > 0){
      big.drops.indx <- max(big.drops)
      tip.ini <- inis[big.drops.indx]
      tip.end <- ends[big.drops.indx]
      #abline(v=c(tip.ini,tip.end),lwd=4)
    } else{
      warning("No big chunks of fluorescence decrease!")
      tip.ini <- tip.end <- which.min(fluo.diff)
    }
  } else{
    warning("No fluorescence decreases on the time slice!")
    tip.ini <- tip.end <- which.min(fluo.diff)
  }
  
  fluo.indxs <- tip.ini:tip.end
  
  if(length(fluo.indxs) > 1){
    sharpest.slope <- fluo.indxs[which.min(diff(t.slice[fluo.indxs]))]
  } else{sharpest.slope <- fluo.indxs}
  #SLOW VERSION - FIND MIN INTERCEPT IN MULTIPLE FITS
  fit.vec <- (sharpest.slope - n.pts):(sharpest.slope + 1)
  
  slope.sharpness <- c()
  for(i in fit.vec){
    slope.sharpness <-c(slope.sharpness, 
                        sum(diff(t.slice[i:(i + n.pts)])))
  }
  
  fit.indxs <- fit.vec[which.min(slope.sharpness)]:(fit.vec[which.min(slope.sharpness)] + (n.pts - 1))
  
  intercept <- mean(fit.indxs - slope * t.slice[fit.indxs])
  return(intercept)
}
#-------------------------------------------------------------------------------
RemoveTipOutliers <- function(tip.loc, dg = 2, spn = 0.2, tol = 3.5){
	tip.loc.coarse <- SmoothWithLoess(tip.loc, spn, dg)
	MZ_score <- abs((0.6745 * ((tip.loc.coarse - tip.loc) - median(tip.loc.coarse - tip.loc))) / mad(tip.loc.coarse - tip.loc))
	if(any(MZ_score > tol)){
		tip.loc[which(MZ_score > tol)] <- tip.loc.coarse[which(MZ_score > tol)]
	}
	return(tip.loc)
}
#-------------------------------------------------------------------------------



