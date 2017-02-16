#-------------------------------------------------------------------------------
FilterWithWavelet <- function(y, smth.vec = c(F, T, T, T, F, F, F)){
  #c(F, T, T, T, F, F, F) is combines lowpass filters by removing higher detail
  #levels with a highpass filter removing noise by not including
  #the lowest level of detail
  lvls <- length(smth.vec) # levels number and their choice controls smoothing
  y.mra <- wavelets::mra(y, filter = "la8", n.levels = lvls, 
                         boundary = "periodic",method = "modwt")
  
  y.osc <- rep(0, length(y))
  for(i in 1:lvls){
    if(smth.vec[i]){y.osc <- y.osc + y.mra@D[[i]]}
  }
  
  trend.ind <- min(c(max(which(smth.vec == T)), length(smth.vec)))
  y.trend <- y.mra@S[[trend.ind]]
  
  smth.ind <- max(min(which(smth.vec == T)) - 1, 1)
  y.smth <- y.mra@S[[smth.ind]]
  
  return(cbind.data.frame("smth" = y.smth, "osc" = y.osc, "trend" = y.trend))
}
#-------------------------------------------------------------------------------
Filter2DWavelet <- function(imaj, smth.vec = c(F, T, T, T, F, F, F)){ #implement plot version
  lvls <- length(smth.vec)
  imaj.mra <- waveslim::mra.2d(imaj, wf = "la8", J = lvls, method = "modwt",
                               boundary = "periodic")
  selected.level <- grepl("HL", names(imaj.mra)) #LLJ for Smooth approx. at level J
  out.imaj <- matrix(0, nrow = dim(imaj)[1], ncol = dim(imaj)[2])
  i <- 1
  j <- 1
  for(detail in imaj.mra){
    if(selected.level[i]){
      if(smth.vec[j]){out.imaj <- out.imaj + detail}
      j <- j + 1
    }
    i <- i + 1
  }
  colnames(out.imaj) <- colnames(imaj)
  rownames(out.imaj) <- rownames(imaj)
  return(out.imaj)
}
#-------------------------------------------------------------------------------