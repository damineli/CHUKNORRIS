#-------------------------------------------------------------------------------
ExtractRidges <- function(wvlt, filter.ridge = "coi", do.summary = FALSE){
  #TODO: RE-implement summary 
  #According to biwavelet plotting function
  if (wvlt$type == "xwt") {
    zvals <- log2(wvlt$power) / (wvlt$d1.sigma * wvlt$d2.sigma)
    get.delay <- TRUE
  } else if (wvlt$type == "wtc" | wvlt$type == "pwtc") {
    zvals <- wvlt$rsq
    wvlt$signif <- wvlt$rsq/wvlt$signif
    get.delay <- TRUE
  } else {
    zvals <- log2(abs(wvlt$power/wvlt$sigma2))
    get.delay <- FALSE
  }
  
  # Single peak corresponding to max power
  ridge.max <- GetMaxPowerRidges(wvlt, pwr = zvals, filter.ridge = filter.ridge, 
                                 get.delay = get.delay, do.summary = do.summary)
  # Multiple power peaks filtered by significance:"signif" coi:"coi" or none:""
  ridge.pwr.pks <- GetPowerPeakRidges(wvlt, pwr = zvals, 
                                      filter.ridge = filter.ridge, 
                                      get.delay = get.delay,
                                      do.summary = do.summary)
  
  ridges <- list("max" = ridge.max, "pwr.pks" = ridge.pwr.pks)
  return(ridges)
}
#-------------------------------------------------------------------------------
GetMaxPowerRidges <- function(wvlt, pwr, p = 1, filter.ridge = "signif", 
                              get.delay = FALSE, do.summary = FALSE){
  
  t.points <- seq(dim(pwr)[2]) 
  y.indxs <- apply(pwr, 2, which.max)
  #plot(wvlt$xaxis,wvlt$period[y.indxs],xlab="Time (min)", ylab="Period (min)")
  
  pwr.rdg <- sapply(1:length(y.indxs), GetWvltValue, tbl = pwr, 
                    ln.indxs = y.indxs) 
  
  sig.rdg <- sapply(1:length(y.indxs), GetWvltValue, tbl = wvlt$signif, 
                    ln.indxs = y.indxs)
  
  if(filter.ridge == "signif"){
    y.indxs[sig.rdg < p] <- NA
    #points(wvlt$xaxis,wvlt$period[y.indxs],col="red",pch=19,cex=0.5)
  } else if(filter.ridge == "coi"){
    y.indxs[sig.rdg < p] <- NA
    in.coi <- wvlt$period[y.indxs] > wvlt$coi[t.points]
    y.indxs[in.coi] <- NA
    #points(wvlt$xaxis,wvlt$period[y.indxs],col="green",pch=19,cex=0.5)
  }
  #points(y=log2(wvlt$period[y.indxs]), x=growth.time)
  #hist(log2(wvlt$period[y.indxs]))
  
  if(!all(is.na(y.indxs))){
    per.rdg <- wvlt$period[y.indxs]
    wave.rdg <- sapply(1:length(y.indxs), GetWvltValue, tbl = wvlt$wave, 
                       ln.indxs = y.indxs) #real part -> amplitude
    phs.rdg <- Arg(wave.rdg)
    amp.rdg <- Re(wave.rdg)
    env.rdg <- Mod(wave.rdg) #actual amplitude envelope
    # (Re(Wave[s.ind, ])/sqrt(Scale[s.ind])) * 
    #   dj * sqrt(dt)/(pi^(-1/4) * 0.776)
    #plot(wave.rdg,type="l")
    #lines(tip.ts$growth.vel.wvlt,col="red",lwd=2)
    #lines(phs.rdg,type="l",col="darkgreen")
    
    ridges.dat <- cbind.data.frame("per" = per.rdg, "phs" = phs.rdg, 
                                   "amp" = amp.rdg, "env" = env.rdg, 
                                   "sig" = sig.rdg, "pwr" = pwr.rdg)
    
    if(get.delay){
      delay.rdg <- per.rdg * (phs.rdg / (2 * pi))
      ridges.dat['delay'] <- delay.rdg
    }
    
    if(do.summary){
      ridges.summary <- SummarizeWvltRidge(ridges.dat)
    } else{ridges.summary <- NULL}
    
  } else{ridges.dat <- ridges.summary <- wave.rdg <- NA} # Change to NULL
  
  return(list("ridges" = ridges.dat, "wave" = wave.rdg, 
              "summary" = ridges.summary))
}
#-------------------------------------------------------------------------------
GetWvltValue <- function(i, tbl, ln.indxs){return(tbl[ln.indxs[i], i])}
#-------------------------------------------------------------------------------
GetPowerPeakRidges<-function(wvlt, pwr, filter.ridge = "signif", 
                             get.delay = FALSE, do.summary = FALSE){
  #filter.ridge.opts = c("signif","coi","")
  
  t.points <- seq(dim(pwr)[2])
  
  if(filter.ridge == "coi"){
    pks.indx.list <- sapply(t.points, FindSigRidgePeaksOutCOI, wvlt = wvlt,
                            pwr = pwr)
  } else if(filter.ridge == "signif"){
    pks.indx.list <- sapply(t.points, FindSigRidgePeaks, wvlt = wvlt, pwr = pwr)
  } else{
    pks.indx.list <- sapply(t.points, FindRidgePeaks, wvlt = wvlt, pwr = pwr)
  }
  
  max.lngth <- max(unlist(lapply(pks.indx.list, length)))
  
  if(max.lngth != 0){
    #Sort by power!
    l <- lapply(pks.indx.list, FillWithNA, max.lngth = max.lngth)
    
    df.pk.indxs <- data.frame(matrix(unlist(l), nrow=(dim(pwr)[2]), byrow=T))
    
    per.tbl <- apply(df.pk.indxs, 2, GetWvltPerLine, per = wvlt$period)
    phs.tbl <- apply(df.pk.indxs, 2, GetWvltPerLine, per = wvlt$phase) #to calculate delays after
    
    wave.tbl <- apply(df.pk.indxs, 2, GetWvltValsLine, tbl = wvlt$wave)
    wave.sum <- apply(wave.tbl, 1, sum, na.rm = T)
    
    phs.rdg <- Arg(wave.sum)
    amp.rdg <- Re(wave.sum)
    env.rdg <- Mod(wave.sum) #actual amplitude
    
    sig.tbl <- apply(df.pk.indxs, 2, GetWvltValsLine, tbl = wvlt$signif)
    pwr.tbl <- apply(df.pk.indxs, 2, GetWvltValsLine, tbl = pwr)
    
    ridges <- list("indxs" = df.pk.indxs, "per" = per.tbl, "phs.vec" = phs.rdg, 
                   "phs" = phs.tbl, "amp" = amp.rdg, "env" = env.rdg,
                   "sig" = sig.tbl, "pwr" = pwr.tbl, "wave" = wave.tbl)
    
    if(get.delay){
      delay.tbl <- per.tbl * (phs.tbl / (2 * pi))
      ridges <- c(ridges, delay = list(delay.tbl))
    }
    
    if(do.summary){
      ridges.summary <- SummarizeWvltRidge(ridges)
    } else{ridges.summary <- NULL}
    
  } else{ridges <- ridges.summary <- wave.sum <- NA} # change for NULL
  
  return(list("ridges" = ridges, "wave" = wave.sum, "summary" = ridges.summary))
}
#-------------------------------------------------------------------------------
FindRidgePeaks <- function(t.point, wvlt, pwr, p=1){
  pk.indxs <- which(diff(sign(diff(pwr[, t.point]))) == -2) + 1
  return(pk.indxs)
}
#-------------------------------------------------------------------------------
FindSigRidgePeaks <- function(t.point, wvlt, pwr, p=1){
  pk.indxs <- which(diff(sign(diff(pwr[, t.point]))) == -2) + 1
  is.signif <- wvlt$signif[pk.indxs, t.point] >= p
  return(pk.indxs[is.signif])
}
#-------------------------------------------------------------------------------
FindSigRidgePeaksOutCOI <- function(t.point, wvlt, pwr, p=1){
  pk.indxs <- which(diff(sign(diff(pwr[, t.point]))) == -2) + 1
  is.signif <- wvlt$signif[pk.indxs, t.point] >= p
  not.in.coi <- wvlt$period[pk.indxs] < wvlt$coi[t.point]
  return(pk.indxs[is.signif & not.in.coi])
}
#-------------------------------------------------------------------------------
FillWithNA <- function(vec, max.lngth){
  return(sort(c(vec, rep(NA, max.lngth - length(vec))), decreasing = T, 
              na.last = T))
}
#-------------------------------------------------------------------------------
GetWvltPerLine <- function(y.indxs, per){return(per[y.indxs])}
#-------------------------------------------------------------------------------
GetWvltValsLine <- function(y.indxs, tbl){return(sapply(1:length(y.indxs), 
                                                        GetWvltValue, tbl = tbl, 
                                                        ln.indxs = y.indxs))}
#-------------------------------------------------------------------------------