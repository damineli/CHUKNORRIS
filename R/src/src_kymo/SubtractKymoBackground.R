#-------------------------------------------------------------------------------
SubtractKymoBackground <- function(kymo, tip.loc, qntl = 0.99, tip.mrgn = 5){
  bg.val <- CalculateBackground(tip.loc, kymo, qntl = qntl, tip.mrgn = tip.mrgn)
  kymo.clean <- kymo - bg.val
  kymo.out <- apply(kymo.clean, 2, FromNegativeToZero)
  return(kymo.out)
}
#-------------------------------------------------------------------------------
CalculateBackground <- function(tip.loc, imaj, qntl = 0.9, tip.mrgn = 5, 
                                fg = FALSE){
  exceed <- which(floor(tip.loc + tip.mrgn) > dim(imaj)[2])
  if(length(exceed) > 0){
    row.limit <- min(which(floor(tip.loc + tip.mrgn) > dim(imaj)[2])) - 1
    message(paste("Kymograph frames", (row.limit + 1), "to", dim(imaj)[2], 
                  "excluded"))
  } else{row.limit <- length(tip.loc)}
  
  if(!fg){
    # Return only the background after a margin
    imaj.nw <- t(sapply(1:row.limit, GetNewFracRowBg, imaj = imaj, 
                        tip.loc = tip.loc + tip.mrgn))
  } else{
    # Return only the foreground after a margin
    imaj.nw <- t(sapply(1:row.limit, GetNewFracRowFg, imaj = imaj, 
                        tip.loc = tip.loc + tip.mrgn))
  }
  #image.plot(imaj.nw)
  bg.val <- quantile(c(imaj.nw), probs = qntl, na.rm = T)
  
  return(bg.val)
}
#-------------------------------------------------------------------------------
GetNewFracRowBg <- function(row.n, imaj, tip.loc){
  tip.round <- floor(tip.loc[row.n])
  tip.frac <- (tip.loc%%1)[row.n]
  #if(length(2:tip.round) != 0){}
  ind.vec <- tip.round:dim(imaj)[2]
  zvals.vec <- imaj[row.n, ]
  new.row <- sapply(ind.vec, CalcNewFracVal, frac = tip.frac, vec = zvals.vec)
  val <- c(rep(NA, (dim(imaj)[2]) - length(new.row)), new.row)
  
  return(val)
}
#-------------------------------------------------------------------------------
GetNewFracRowFg <- function(row.n, imaj, tip.loc){
  tip.round <- floor(tip.loc[row.n])
  tip.frac <- (tip.loc%%1)[row.n]
  #if(length(2:tip.round) != 0){}
  ind.vec <- 2:tip.round
  zvals.vec <- imaj[row.n, ]
  new.row <- sapply(ind.vec, CalcNewFracVal, frac = tip.frac, vec = zvals.vec)
  val <- c(rep(NA, (dim(imaj)[2]) - length(new.row)), new.row)
  
  return(val)
}
#-------------------------------------------------------------------------------
FromNegativeToZero <- function(vec){
  vec[which(vec < 0)] <- 0
  return(vec)
}
#-------------------------------------------------------------------------------

