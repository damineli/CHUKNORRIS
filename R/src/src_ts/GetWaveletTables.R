#-------------------------------------------------------------------------------
GetWaveletTable <- function(dat.all, max.ridges, pwr.pks, fl.nm, var.nm, #change full.id to fl.nm
                            max.comp = 5){
  if(is.null(dim(max.ridges))){
    max.ridges <- data.frame(matrix(NA, nrow = dim(dat.all)[1], ncol = 6))
  }
  rownames(dat.all) <- row.names(max.ridges)
  tbl.1 <- cbind(dat.all, max.ridges)
  tbl.1.nms <- c("time", "raw", "smth", "osc", "trend", "per", "phs", "amp", "env", 
                 "sig", "pwr")
  
  null.col <- rep(NA, dim(pwr.pks$per)[1])
  prev.lines <- null.col
  pwr.pks.col.nms <- c()
  
  for(i in 1:max.comp){
    if(i <=  dim(pwr.pks$per)[2]){
      current.line <- cbind(pwr.pks$per[, i], pwr.pks$phs[, i], pwr.pks$sig[, i],
                            pwr.pks$pwr[, i])
    } else{
      current.line <- cbind(null.col, null.col, null.col, null.col)
    }
    prev.lines <- cbind(prev.lines, current.line)
    pwr.pks.col.nms <- c(pwr.pks.col.nms, paste(c("per", "phs", "sig", "pwr"),
                                                ".", i, sep=""))
  }
  
  wvlt.tbl <- cbind(tbl.1, prev.lines[, 2:dim(prev.lines)[2]], pwr.pks$phs.vec,
                    pwr.pks$amp, pwr.pks$env)
  all.col.names <- c(tbl.1.nms, pwr.pks.col.nms, "phs.tot", "amp.tot", "env.tot")
  colnames(wvlt.tbl) <- all.col.names
  
  wvlt.tbl <- AddNameColumn(wvlt.tbl, var.nm)
  wvlt.tbl <- AddNameColumn(wvlt.tbl, fl.nm) #change full.id to fl.nm
  colnames(wvlt.tbl)[1:2] <- c("fl.nm", "var.nm") #change id to fl.nm
  
  return(wvlt.tbl)
}
#-------------------------------------------------------------------------------
GetJointWaveletTable <- function(ts1, ts2, var1.nm = "var1", var2.nm = "var2",
                                 fl.nm = "", technique = "?", max.ridges, 
                                 pwr.pks, max.comp = 5){
  
  dat.all <- cbind(ts1, ts2[, 2])
  if(is.null(dim(max.ridges))){
    max.ridges <- data.frame(matrix(NA, nrow = dim(dat.all)[1], ncol = 6))
  }
  rownames(dat.all) <- row.names(max.ridges)
  tbl.1 <- cbind(dat.all, max.ridges)
  tbl.1.nms <- c("time", "var1", "var2", "per", "phs", "amp", "env", "sig", "pwr", "delay")
  
  null.col <- rep(NA, dim(pwr.pks$per)[1])
  prev.lines <- null.col
  pwr.pks.col.nms <- c()
  
  for(i in 1:max.comp){
    if(i <=  dim(pwr.pks$per)[2]){
      current.line <- cbind(pwr.pks$per[, i], pwr.pks$phs[, i], pwr.pks$sig[, i], pwr.pks$pwr[, i], pwr.pks$delay[, i])
      #print(dim(current.line))
    } else{
      current.line <- cbind(null.col, null.col, null.col, null.col, null.col)
    }
    prev.lines <- cbind(prev.lines, current.line)
    pwr.pks.col.nms <- c(pwr.pks.col.nms, paste(c("per", "phs", "sig", "pwr", "delay"),".", i,sep=""))
  }
  
  wvlt.tbl <- cbind(tbl.1, prev.lines[, 2:dim(prev.lines)[2]], pwr.pks$phs.vec, pwr.pks$amp, pwr.pks$env)
  all.col.names <- c(tbl.1.nms, pwr.pks.col.nms, "phs.tot", "amp.tot", "env.tot")
  colnames(wvlt.tbl) <- all.col.names
  
  wvlt.tbl <- AddNameColumn(wvlt.tbl, technique)
  wvlt.tbl <- AddNameColumn(wvlt.tbl, var2.nm)
  wvlt.tbl <- AddNameColumn(wvlt.tbl, var1.nm)
  wvlt.tbl <- AddNameColumn(wvlt.tbl, fl.nm)
  colnames(wvlt.tbl)[1:4] <- c("fl.nm", "var1.nm", "var2.nm", "technique")
  
  return(wvlt.tbl)
}

#-------------------------------------------------------------------------------