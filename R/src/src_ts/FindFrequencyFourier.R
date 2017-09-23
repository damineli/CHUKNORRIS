#-------------------------------------------------------------------------------
#----------------------------------------------------FOURIER-BASED METHODS BLOCK
#-------------------------------------------------------------------------------
FindFrequencyFourier <- function(ts, sf){
  ##############################################################################
  ## TODO: Calculate CI, get amplitude and phase (as wvlt?) extend to coherence
  ##############################################################################
  min.f <- 1 / ((length(ts) * (1 / sf)) / 2)
  frq.spec.pgram <- spec.pgram(ts, log="no", plot = FALSE, fast = FALSE)
  frq.spec.ar <- spec.ar(ts, log="no", plot = FALSE)  
  spec.pks <- list("pgram" = FindFreqPeak(frq.spec.pgram, sf = sf, min.f),
                   "ar" = FindFreqPeak(frq.spec.ar, sf = sf, min.f))
  
  spec.info <- list("spec.pgram" = frq.spec.pgram,
                    "spec.ar" = frq.spec.pgram, "spec.pks" = spec.pks)
  
  return(spec.info)
}
#-------------------------------------------------------------------------------
FindFreqPeak <- function(frq.spec, sf, min.f = 0){
  ##############################################################################
  ## TODO: Better peak finding algorithm (interpolation)
  ##############################################################################
  spec.frq <- frq.spec$freq * sf
  spec.pwr <- frq.spec$spec
  spec.pwr[spec.frq < min.f] <- NA
  
  frq.pk <- spec.frq[which.max(spec.pwr)]
  if(length(frq.pk) == 0){frq.pk <- NA}
  #abline(h=log2(1/(tst$freq[which.max(tst$spec)]*(1/(4/60)))))
  #per.pk <- (1 / frq.pk) / 60 # Period in minutes
  #spec.pks = c("frq" = frq.pk, "per" = per.pk)
  #lsp(x=growth.dat$raw.osc[-1],times=time.vec[-1],from=0.2,to=2,type="period",alpha = 0.05)->lmb
  return(frq.pk)
}
#-------------------------------------------------------------------------------
