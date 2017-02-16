#-------------------------------------------------------------------------------
AddIfNotNull <- function(tbl.1, tbl.2){
  if((is.null(tbl.1)) & (is.null(tbl.2))){
    return(NULL)
  } else if(is.null(tbl.1)){
    return(tbl.2)
  } else if(is.null(tbl.2)){
    return(tbl.1)
  } else if((!is.null(tbl.1)) & (!is.null(tbl.2))){
    return(rbind(tbl.1, tbl.2))
  } else{message("Conflicting table content!")}
}
#-------------------------------------------------------------------------------
AddNameColumn <- function(dat, nm){return(cbind(rep(nm, dim(dat)[1]), dat))}
#-------------------------------------------------------------------------------
