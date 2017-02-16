#-------------------------------------------------------------------------------
GetFluoTimeSeries <- function(kymo, ini.ind, avg.width){
  return(apply(kymo[, ini.ind:(ini.ind + avg.width)], 1, median))
}
#-------------------------------------------------------------------------------