#-------------------------------------------------------------------------------
GetNewFracRow <- function(row.n, imaj, tip.loc){
  tip.round <- floor(tip.loc[row.n])
  tip.frac <- (tip.loc%%1)[row.n]
  #if(length(2:tip.round) != 0){}
  ind.vec <- 2:tip.round
  zvals.vec <- imaj[row.n, ]
  new.row <- sapply(ind.vec, CalcNewFracVal, frac = tip.frac, vec = zvals.vec)
  val <- c(rep(0, (dim(imaj)[2]) - length(new.row)), new.row)
  
  return(val)
}
#-------------------------------------------------------------------------------
CalcNewFracVal <- function(i, frac, vec){
  return(weighted.mean(c(vec[i], vec[i - 1]), c(frac, 1 - frac)))
}
#-------------------------------------------------------------------------------
