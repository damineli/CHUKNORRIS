#-------------------------------------------------------------------------------
DoContWvlt <- function(ts, dj = 1 / 100, wvlt.type = "morlet", pad = F, 
                       p = 0.95, filter.ridge = "coi", do.summary = FALSE){
  ts.cwt <- biwavelet::wt(ts, pad = pad, dj = dj, mother = wvlt.type, 
                          sig.level = p)
  ridges <- ExtractRidges(ts.cwt, filter.ridge, do.summary)
  
  return(list("cwt" = ts.cwt,
              "ridges" = ridges))
}
#-------------------------------------------------------------------------------
DoCrossWvlt <- function(t1, t2, dj = 1 / 100, wvlt.type = "morlet", pad = F, 
                        p = 0.95, filter.ridge = "coi", do.summary = FALSE){
  ts.xwt <- biwavelet::xwt(t1, t2, pad = pad, dj = dj, mother = wvlt.type, 
                          sig.level = p)
  ridges <- ExtractRidges(ts.xwt, filter.ridge, do.summary)
  return(list("wvlt" = ts.xwt,
              "ridges" = ridges))
}
#-------------------------------------------------------------------------------
DoWvltCoher <- function(t1, t2, dj = 1 / 100, wvlt.type = "morlet", pad = F, 
                        p = 0.95, filter.ridge = "coi", do.summary = FALSE, 
                        nrands = 300){
  wvlt.coher <- biwavelet::wtc(t1, t2, pad = pad, dj = dj, mother = wvlt.type, 
                               sig.level = p, nrands = nrands)
  ridges <- ExtractRidges(wvlt.coher, filter.ridge, do.summary)
  return(list("wvlt" = wvlt.coher,
              "ridges" = ridges))
  }
#-------------------------------------------------------------------------------