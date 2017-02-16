#-------------------------------------------------------------------------------
AlignByTip <- function(tip.loc, imaj){
  exceed <- which(tip.loc > dim(imaj)[2])
  if(length(exceed) > 0){
    row.limit <- min(which(tip.loc > dim(imaj)[2])) - 1
    message(paste("Kymograph frames", (row.limit + 1), "to", dim(imaj)[2], 
                  "excluded"))
  } else{row.limit <- length(tip.loc)}
  
  imaj.nw <- t(sapply(1:row.limit, GetNewFracRow, imaj = imaj,
                      tip.loc = tip.loc))
  
  return(imaj.nw[, dim(imaj.nw)[2]:1])
}
#-------------------------------------------------------------------------------