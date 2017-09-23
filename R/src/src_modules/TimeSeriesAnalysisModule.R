#-------------------------------------------------------------------------------
AnalyzeTimeSeries <- function(dat, ts.par, filter.par, wvlt.par, out.par, 
                              var.nm = "x", fl.nm = "file.x"){
  # Create holes if NA
  dat <- dat[(!is.na(dat[, 1])) & (!is.na(dat[, 2])), ]
  # Open plotting device
  if(out.par$do.plot){
    pdf(paste(out.par$out.path, "OscillationAnalysis_", fl.nm, ".pdf", sep = ""), 
        width = 7.5, height = 6)
    par(mar = c(5, 4.5, 4, 2) + 0.1)
  }
  
  dat.raw <- dat # save raw data to output later
  if(ts.par$time.unit == "min"){
    # Convert to seconds (for more precise calculations)
    dat[, 1] <- dat[, 1] * 60
    filter.par$low.per <- filter.par$low.per * 60
    filter.par$high.per <- filter.par$high.per * 60
    }
  # Pre-filter if selected
  if(ts.par$do.prefilter){dat <- PreFilterTimeSeries(dat, ts.par, out.par)}
  # Get time step
  time.step <- GetSamplingIntervalMode(dat[, 1], dg = ts.par$time.dgts)
  # Calculate sampling frequency
  sampling.freq <- 1 / time.step
  # If frequency bands where not selected, calculate given specified intervals
  if(is.null(filter.par$smth.vec)){
    filter.par$smth.vec <- GetSmoothVec(sampling.freq,
                                        series.length = dim(dat)[1] * time.step,
                                        low.per = filter.par$low.per, 
                                        high.per = filter.par$high.per)
  }
  
  
  # Filter with loess first
  if(is.null(filter.par$loess.first)){filter.par$loess.first <- FALSE}
  if(filter.par$loess.first){
    trnd <- SmoothWithLoess(dat[, 2], spn = filter.par$loess.spn, 
                            degree = filter.par$loess.degree)
    dat.filt <- cbind(dat, 
                      FilterWithWavelet(dat[, 2] - trnd, 
                                        smth.vec = filter.par$smth.vec))
    colnames(dat.filt)[1:2] <- c("time", "raw")
    dat.filt$trend <- dat.filt$trend + trnd
    dat.filt$smth <- dat.filt$smth + trnd
    
  } else{
    # Filter with Wavelet or skip
    if(length(filter.par$smth.vec) == 1){
      if(filter.par$smth.vec == "skip"){
        trnd <- SmoothWithLoess(dat[, 2], spn = 0.4, degree = 1)
        dat.filt <- cbind(dat, cbind.data.frame("smth" = dat[, 2], 
                                                "osc" = dat[, 2] - trnd, 
                                                "trend" = trnd))
        colnames(dat.filt)[1:2] <- c("time", "raw")
      }
    } else{
      # Filter with a discrete wavelet transform and join all time series data  
      dat.filt <- cbind(dat, 
                        FilterWithWavelet(dat[, 2], smth.vec = filter.par$smth.vec))
      colnames(dat.filt)[1:2] <- c("time", "raw")
    }
    
  }
   
  # Get table with all summary statistics and Fourier transform peak
  summary.stats <- GetAllSummaryStats(dat.filt, var.nm, fl.nm)
  if(ts.par$time.unit == "min"){
    # Convert back to minutes
    summary.stats$period.pgram <- summary.stats$period.pgram / 60
    summary.stats$period.ar <- summary.stats$period.ar / 60
    dat.filt[, 1] <- dat.filt[, 1] / 60
    
  }
  # Analyze oscillations with continuous wavelet and get results table
  osc.analysis <- AnalyzeOscillations(dat.filt, wvlt.par, var.nm, fl.nm)
  
  # Temporary plots
  # Plot filtering results
  if(out.par$do.plot){
    plot(dat.filt[ ,1], dat.filt$raw, type="o", col = "grey", lwd = 1.5, 
         xlab = paste("Time (", ts.par$time.unit, ")", sep = ""), 
         ylab = out.par$ylab, main = "Filtering results")
    lines(dat.filt[ ,1], dat.filt$trend, 
          col = adjustcolor("black", alpha.f = 1))
    lines(dat.filt[ ,1], dat.filt$smth, 
          col = adjustcolor("red", alpha.f = 0.75))
    legend("bottomleft", legend = c("Prefiltered", "Filtered", "Trend"), 
           col = c("grey", "red","black"), bty = "n", lwd = 1, 
           pch = c(1, -1 , -1))
    
    plot(dat.filt[ ,1], dat.filt$osc, type="l", col = "blue", lwd = 2,
         xlab = paste("Time (", ts.par$time.unit, ")", sep = ""), 
         ylab = out.par$ylab, main = "Filtered and detrended series")
    
  }
  
  # Plot wavelet
  if(out.par$do.plot){
    wvlt.pal <- colorRampPalette(tim.colors(256), bias = 0.6)
    plot(osc.analysis$cwt, type = "power.corr.norm", 
         main = "Continuous Wavelet Transform",  
         xlab = paste("Time (", ts.par$time.unit, ")", sep = ""),  
         ylab = paste("Period (", ts.par$time.unit, ")", sep = ""),
         fill.cols = wvlt.pal(512))
    for(i in 1:dim(osc.analysis$ridges$pwr.pks$ridges$per)[2]){
      points(dat.filt[, 1], log2(osc.analysis$ridges$pwr.pks$ridges$per[, i]), 
             col = "white", pch = 19, cex = 0.5)
    }
    abline(h = log2((summary.stats$period.pgram[3])), col = "white",
           lty = "dotdash", lwd = 2)
    legend("bottomleft", c("Power peaks (ridges)", "Fourier peak"),
           pch = 19, pt.cex = c(0.5, 0), col = "white", lty = "dotdash", 
           lwd = c(-1,2), text.col = "white")
  }
  
 if(wvlt.par$do.summary){
 	  # GMM
 	  if(!is.null(osc.analysis$ridges$max$ridges) & (length(osc.analysis$ridges$max$ridges) != 1)){
 	  	gmm.max <- SummarizeWvltRidge(osc.analysis$ridges$max$ridges)
 	  	gmm.pwr.pks <- SummarizeWvltRidge(osc.analysis$ridges$pwr.pks$ridges)
 	  	if(!is.null(gmm.max$per$gmm$mu)){
 	  		PlotGMM(gmm.max$per$gmm, obs = c(osc.analysis$ridges$max$ridges$per), 
 	  		xlb = "Period (min)", ttl = paste(var.nm, "max ridge"), cte = 1.75)
 	  		}
 	  		if(!is.null(gmm.pwr.pks$per$gmm$mu)){
 	  			PlotGMM(gmm.pwr.pks$per$gmm, obs = c(osc.analysis$ridges$pwr.pks$ridges$per), 
 	  			xlb = "Period (min)", ttl = paste(var.nm, "all ridges"), cte = 1.75)
 	  			}
 	  		}
 	  	}

  # Close plotting device
  if(out.par$do.plot){
    dev.off()
  }

  out.list <- list("raw" = dat.raw,
                   "cwt" = osc.analysis$cwt,
                   "wvlt.tbl" = osc.analysis$wvlt.tbl,
                   "summary.stats.tbl" = summary.stats)
  
  # Save analysis tables
  if(out.par$do.save.data){
    write.csv(out.list$wvlt.tbl, 
              paste(out.par$out.path, "OscillationAnalysisTbl_", fl.nm, ".csv", 
                    sep = ""), row.names = FALSE)
    write.csv(out.list$summary.stats.tbl, 
              paste(out.par$out.path, "SummaryStats_", fl.nm, ".csv", sep = ""), 
              row.names = FALSE)
  }
  
  return(out.list)
}
#-------------------------------------------------------------------------------
AnalyzeOscillations <- function(dat, wvlt.par, var.nm, fl.nm){
  cwt.lst <- DoContWvlt(cbind(dat$time, dat$osc),
                       dj = wvlt.par$dj,
                       wvlt.type = wvlt.par$wvlt.type,
                       pad = wvlt.par$do.pad,
                       p = wvlt.par$p,
                       filter.ridge = wvlt.par$filter.ridge,
                       do.summary = wvlt.par$do.summary)
  ts.cwt <- cwt.lst$cwt
  ridges <- cwt.lst$ridges
  
  wvlt.tbl <- GetWaveletTable(dat.all = dat,
                              max.ridges = ridges$max$ridges,
                              pwr.pks = ridges$pwr.pks$ridges,
                              fl.nm = fl.nm,
                              var.nm = var.nm,
                              max.comp = wvlt.par$max.comp)
  
  return(list("cwt" = ts.cwt,
              "wvlt.tbl" = wvlt.tbl,
              "wvlt.summary" = ridges$summary,
              "ridges" = ridges))
}
#-------------------------------------------------------------------------------