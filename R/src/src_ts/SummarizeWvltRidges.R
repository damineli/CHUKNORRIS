#-------------------------------------------------------------------------------
SummarizeWvltRidge <- function(ridges){
  
  per.vec <- as.vector(unlist(ridges$per)) # since values are pre-filtered...
  is.sig <- (!is.na(per.vec))
  per.vec <- per.vec[is.sig]
  
  sig.vec <- as.vector(unlist(ridges$sig))[is.sig]
  wts <- (sig.vec / sum(sig.vec))
  phs.vec <- as.vector(unlist(ridges$phs))[is.sig]
  env.vec <- ridges$env[is.sig[1:length(ridges$env)]]
  wts.env <- sig.vec[1:length(env.vec)] / sum(sig.vec[1:length(env.vec)])
  
  per.dstr <- FitDistr(obs = per.vec, wts = wts)
  phs.dstr <- FitDistr(obs = phs.vec, wts = wts, fit.gmm = FALSE)
  amp.dstr <- FitDistr(obs = env.vec, wts = wts.env)
  
  wvlt.rdg.dstrbs <- list("per" = per.dstr, "phs" = phs.dstr, "amp" = amp.dstr)
  
  if("delay"%in%names(ridges)){
    delay.vec <- as.vector(unlist(ridges$delay))[is.sig]
    delay.dstr <- FitDistr(obs = delay.vec, wts = wts)
    wvlt.rdg.dstrbs <- c(wvlt.rdg.dstrbs , delay = list(delay.dstr))
  }
  
  return(wvlt.rdg.dstrbs)
}
#-------------------------------------------------------------------------------
FitDistr <- function(obs, wts, adjst = 1, max.pks = 4, min.pk.d = 0.01, 
                     min.obs = 5, fit.gmm = TRUE){
  if(length(obs) > min.obs){
    # Fit a weighted density curve
    obs.d <- density(obs, adjust = adjst, weights = wts)
    #d <- density(obs, adjust = adjst)
    
    # Find density peaks
    d.pks <- FindDensPeaks(d = obs.d)
    min.d <- (max(obs.d$y[match(d.pks, obs.d$x)])*min.pk.d)
    d.pks <- d.pks[(obs.d$y[match(d.pks, obs.d$x)] >= min.d)]
    max.pks <- min(max.pks, floor(length(obs)/min.obs))
    d.pks.gmm <- d.pks[1:min(c(max.pks, length(d.pks)))]
    # Require package 'mixtools'
    # Later sort mus by lambda!
    if(fit.gmm & (!any(is.na(d.pks.gmm)))){
      obs.gmm <- RobustGMM(obs, mu = d.pks.gmm)
      obs.gmm.auto <- RobustGMM(obs)
    } else{obs.gmm <- obs.gmm.auto <- NULL}
    
  } else{obs.d <- d.pks <- obs.gmm <- obs.gmm.auto <- NULL}
  
  fit.dstr.lst <- list("density" = obs.d, "dens_pks" = d.pks,
                       "gmm" = obs.gmm, "gmm_auto" = obs.gmm.auto)
  
}
#-------------------------------------------------------------------------------
FindDensPeaks <- function(d){
  peaks.ind <- FindPeaks(d$y)
  max.peaks.ind <- sort(d$y[peaks.ind], decreasing = T, index.return = T)$ix
  #if(length(peaks_ind)>1){d$x[peaks_ind[ max_peaks_ind[1:2]]]}
  #TODO: exclude close peaks?
  return(d$x[peaks.ind[max.peaks.ind]])
}
#-------------------------------------------------------------------------------
FindPeaks <- function(vec){return(which(diff(sign(diff(vec))) == -2) + 1)}
#-------------------------------------------------------------------------------
EstimateDelay <- function(phs, per, wts){
  
  delay <- per * (phs / (2 * pi))
  
  obs.d <- density(obs, adjust = adjst, weights = wts)
  #d <- density(obs, adjust = adjst)
  
  # Find density peaks
  d.pks <- FindDensPeaks(d = obs.d)
  
}
#-------------------------------------------------------------------------------
