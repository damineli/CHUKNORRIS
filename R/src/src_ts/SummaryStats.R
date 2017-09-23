#-------------------------------------------------------------------------------
GetAllSummaryStats <- function(dat, var.nm = "", fl.nm = "", dg = 6){
  # Get time step
  time.step <- GetSamplingIntervalMode(dat[, 1], dg)
  # Calculate sampling frequency
  sampling.freq <- 1 / time.step
  # Analyze all time series
  ts.nms <- c("raw", "smth", "osc", "trend")
  summary.tbl <- c()
  
  for(ts.nm in ts.nms){
    summary.stats <- DoSummaryStats(dat[, ts.nm])
    spec.pks <- FindFrequencyFourier(dat[, ts.nm], 
                                     sf = sampling.freq)$spec.pks
    tbl <- GetStatsSummaryTable(summary.stats, spec.pks, ts.nm = ts.nm, 
                                var.nm = var.nm, fl.nm = fl.nm)
    summary.tbl <- rbind(summary.tbl, tbl) 
  }
  
  return(summary.tbl)
}
#-------------------------------------------------------------------------------
DoSummaryStats <- function(ts.dat){
  
  summary.stats <- list("summary" = summary(ts.dat), 
                        "sd" = sd(ts.dat, na.rm = T), 
                        "mad" = mad(ts.dat, na.rm = T),
                        "iqr.1" = IQR(ts.dat, na.rm = T),
                        "qntls.2" = quantile(ts.dat, na.rm = T, 
                                             prob = c(0.01, 0.99)),
                        "iqr.2" = diff(quantile(ts.dat, na.rm = T, 
                                                prob = c(0.01, 0.99))))
  
  return(summary.stats)
  
}
#-------------------------------------------------------------------------------
GetStatsSummaryTable <- function(summary.stats, spec.pks, ts.nm = "", 
                                 var.nm = "", fl.nm = ""){
  summary.stats <- unlist(summary.stats)
  if(length(summary.stats) == 13){
    summary.stats <- summary.stats[-7] 
  }
  
  stat.dat <- as.data.frame(c(summary.stats, 1 / spec.pks$pgram,
                              1 / spec.pks$ar)) # in seconds by default
  stats.summary.tbl <- t(stat.dat)
  rownames(stats.summary.tbl) <- NULL
  
  s.stats.colnames <- c("min", "25%.quant", "median", "mean", "75%.quant", 
                        "max", "sd", "mad", "iqr.1", "1%.quant", "99%.quant",
                        "iqr.2", "period.pgram", "period.ar")
  
  colnames(stats.summary.tbl) <- s.stats.colnames
  out.tbl <- cbind.data.frame(rep(ts.nm, dim(stats.summary.tbl)[1]), 
                              stats.summary.tbl)
  colnames(out.tbl)[1] <- "ts.nm" 
  
  out.tbl <- cbind.data.frame(rep(var.nm, dim(out.tbl)[1]), out.tbl)
  colnames(out.tbl)[1] <- "var.nm" 
  
  out.tbl <- cbind.data.frame(rep(fl.nm, dim(out.tbl)[1]), out.tbl)
  colnames(out.tbl)[1] <- "fl.nm" 
  
  return(out.tbl)
}
#-------------------------------------------------------------------------------