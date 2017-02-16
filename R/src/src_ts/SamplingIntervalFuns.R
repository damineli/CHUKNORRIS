#-------------------------------------------------------------------------------
GetSamplingIntervalMode <- function(t.vec, dg = 4){
  return(Mode(round(diff(t.vec), digits = dg)))
}
#-------------------------------------------------------------------------------
TestSamplingIntervalHomogeneity <- function(t.vec, dg = 6){
  return(round(sd(round(diff(t.vec), digits = dg)), digits = dg) == 0)
}
#-------------------------------------------------------------------------------
