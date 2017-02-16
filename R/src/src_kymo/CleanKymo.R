#-------------------------------------------------------------------------------
CleanKymo <- function(dirty.kymo){
  # Changes NA/NaN to Zero and Inf to Maximum values of a ratiometric kymograph
  
  # Vectorize matrix
  val.vec <- c(dirty.kymo)
  
  # Compute maximum for non-Na/NaN and finite numbers
  max.val <- max(val.vec[(!is.na(val.vec)) & is.finite(val.vec)])
  
  # Substitutes Na/NaN for 0 and infinite numbers for maximum value
  clean.kymo <- apply(dirty.kymo, 2, SubstituteNaNsAndInfs, max.val)
  
  return(clean.kymo)
}
#-------------------------------------------------------------------------------
SubstituteNaNsAndInfs <- function(vec, max.val){
  # Assigns 0 to NA/NaN and max.val to Inf
  vec[is.na(vec)] <- 0
  vec[!is.finite(vec)] <- max.val
  return(vec)
}
#-------------------------------------------------------------------------------