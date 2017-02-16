#-------------------------------------------------------------------------------
AnalyzeTimeSeriesSynchrony <- function(dat1, dat2, par.lst, wvlt.pair.par, 
                                       out.par, var1.nm = "var1", 
                                       var2.nm = "var2", fl1.nm = "", 
                                       fl2.nm = ""){
  # Analyze series individually
  t1.lst <- AnalyzeTimeSeries(dat = dat1, ts.par = par.lst$par1$ts.par, 
                              filter.par = par.lst$par1$filter.par,
                              wvlt.par = par.lst$par1$wvlt.par, 
                              out.par = par.lst$par1$out.par, 
                              var.nm = var1.nm, fl.nm = fl1.nm)
  
  t2.lst <- AnalyzeTimeSeries(dat = dat2, ts.par = par.lst$par2$ts.par, 
                              filter.par = par.lst$par2$filter.par,
                              wvlt.par = par.lst$par2$wvlt.par, 
                              out.par = par.lst$par2$out.par, 
                              var.nm = var2.nm, fl.nm = fl2.nm)
  
  # Match series pair
  series.pair <- MatchSeries(t1 = cbind(t1.lst$wvlt.tbl$time, 
                                        t1.lst$wvlt.tbl$osc),
                             t2 = cbind(t2.lst$wvlt.tbl$time, 
                                        t2.lst$wvlt.tbl$osc), 
                             smth.n = wvlt.pair.par$interpol.n)
  
  ts1 <- series.pair[[1]]
  ts2 <- series.pair[[2]]
  
  # Define file name for joint analysis
  fl.nm <- paste(fl1.nm, fl2.nm, sep = "_VS_")
  
  # Analyze pair synchornization with Cross Wavelet
  cross.wvlt.lst <- AnalyzePairSynchronization(ts1, ts2, 
                                               var1.nm = var1.nm, 
                                               var2.nm = var2.nm, 
                                               fl.nm = fl.nm, 
                                               wvlt.par = wvlt.pair.par, 
                                               do.wvlt.coher = FALSE)
  # Plotting
  if(out.par$do.plot){
    pdf(paste(out.par$out.path, "Sync_", fl.nm, ".pdf", sep = ""), 
        width = 9, height = 6)
    par(mar = c(5, 4.5, 4, 4.5) + 0.1)
    
    # Plot matched series
    x.rng <- range(c(t1.lst$wvlt.tbl$time, t2.lst$wvlt.tbl$time))
    plot(t1.lst$wvlt.tbl$time, t1.lst$wvlt.tbl$osc, xlim = x.rng, type = "l", 
         col = "slateblue4", lwd = 1, xaxt="n",ylab="",yaxt="n", 
         xlab = paste("Period (", par.lst$par1$ts.par$time.unit, ")"),
         main = "Filtered and matched series")
    lines(ts1, col = adjustcolor("slateblue2", alpha.f = 0.5), lwd = 5)
    axis(2,col.axis = "slateblue4",col = "slateblue4")
    axis(1)
    mtext(par.lst$par1$out.par$ylab, side = 2, line = 3, col = "slateblue4")
    
    par(new = TRUE)
    
    plot(t2.lst$wvlt.tbl$time, t2.lst$wvlt.tbl$osc, xlim = x.rng, type = "l", col = "firebrick4", xaxt="n",ylab="",yaxt="n", 
    xlab = "")
    lines(ts2, col = adjustcolor("firebrick2", alpha.f = 0.5), lwd = 5)
    axis(4,col.axis = "firebrick4",col = "firebrick4")
    mtext(par.lst$par2$out.par$ylab, side = 4, line = 3, col = "firebrick4")
    
    legend("topleft",c(paste("Filtered", var1.nm), paste("Matched", var1.nm), 
                       paste("Filtered", var2.nm), paste("Matched", var2.nm)),
           lwd = c(1, 5, 1, 5), bty = "n",
           col = c("slateblue4", "slateblue2", "firebrick4", "firebrick2"))
    
    # Plot Cross Wavelet
    plot(cross.wvlt.lst$analysis$wvlt, type = "power.norm", 
         plot.phase = TRUE, main = "Cross Wavelet Transform",  
         xlab = paste("Time (", par.lst$par1$ts.par$time.unit, ")", sep = ""),  
         ylab = paste("Period (", par.lst$par1$ts.par$time.unit, ")", sep = ""))
    legend("bottomleft", paste(var1.nm, var2.nm, sep = " vs "), 
           text.col = "white", bg = "transparent")
    }
  if(wvlt.pair.par$do.wvlt.coher){
    # Analyze pair synchornization with Wavelet Coherence
    wvlt.coher.lst <- AnalyzePairSynchronization(ts1 = series.pair[[1]],
                                                 ts2 = series.pair[[2]], 
                                                 var1.nm = var1.nm, 
                                                 var2.nm = var2.nm, 
                                                 fl.nm = fl.nm, 
                                                 wvlt.par = wvlt.pair.par,
                                                 do.wvlt.coher = TRUE)
    if(out.par$do.plot){
      # Plot Wavelet Coherence
      plot(wvlt.coher.lst$analysis$wvlt, type = "power.norm", 
           plot.phase = TRUE, main = "Wavelet Coherence",  
           xlab = paste("Time (", par.lst$par1$ts.par$time.unit, ")", sep = ""),  
           ylab = paste("Period (", par.lst$par1$ts.par$time.unit, ")", sep = ""))
      legend("bottomleft", paste(var1.nm, var2.nm, sep = " vs "), 
             text.col = "white", bg = "transparent")
      }
    
    } else{wvlt.coher.lst <- list("analysis" = NULL, "tbl" = NULL)}
  
  sync.tbl <- AddIfNotNull(cross.wvlt.lst$tbl, wvlt.coher.lst$tbl)
  
  if(out.par$do.plot){dev.off()}
  
  pair.analyses <- list("t1.lst" = t1.lst, "t2.lst" = t2.lst, 
                        "matched.series" = series.pair,
                        "sync.tbl" = sync.tbl,
                        "cross.wvlt" = cross.wvlt.lst$analysis,
                        "wvlt.coher" = wvlt.coher.lst$analysis)
  
  if(out.par$do.save.data){
    write.csv(sync.tbl, paste(out.par$out.path, "SyncTbl_", fl.nm, 
                              ".csv", sep = ""), row.names = FALSE)
  }
  
  return(pair.analyses)
}
#-------------------------------------------------------------------------------
AnalyzePairSynchronization <- function(ts1, ts2, var1.nm = "x", var2.nm = "y", 
                                       fl.nm = "", wvlt.par, 
                                       do.wvlt.coher = FALSE){
  if(do.wvlt.coher){
    analysis <- DoWvltCoher(ts1, ts2, 
                            dj = wvlt.par$dj,
                            wvlt.type = wvlt.par$wvlt.type,
                            pad = wvlt.par$do.pad,
                            p = wvlt.par$p,
                            filter.ridge = wvlt.par$filter.ridge,
                            do.summary = wvlt.par$do.summary,
                            nrand = wvlt.par$nrand)
    tbl <- GetJointWaveletTable(ts1, ts2,
                               var1.nm = var1.nm, var2.nm = var2.nm,
                               fl.nm = fl.nm, technique = "wvlt.coher",
                               max.ridges = analysis$ridges$max$ridges,
                               pwr.pks = analysis$ridges$pwr.pks$ridges,
                               max.comp = wvlt.par$max.comp)
    
  } else{
    analysis <- DoCrossWvlt(ts1, ts2, 
                            dj = wvlt.par$dj,
                            wvlt.type = wvlt.par$wvlt.type,
                            pad = wvlt.par$do.pad,
                            p = wvlt.par$p,
                            filter.ridge = wvlt.par$filter.ridge,
                            do.summary = wvlt.par$do.summary)
    tbl <- GetJointWaveletTable(ts1, ts2,
                               var1.nm = var1.nm, var2.nm = var2.nm,
                               fl.nm = fl.nm, technique = "cross.wvlt",
                               max.ridges = analysis$ridges$max$ridges,
                               pwr.pks = analysis$ridges$pwr.pks$ridges,
                               max.comp = wvlt.par$max.comp)
  }
  return(list("analysis" = analysis, "tbl" = tbl))
}
#-------------------------------------------------------------------------------