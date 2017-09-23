#-------------------------------------------------------------------------------
AnalyzeRatiometricKymograph <- function(k1, k2, tip.find.par, image.par, ts.par,
                                        filter.par, wvlt.par, out.par, var1.nm, 
                                        var2.nm, fl.nm = ""){
  
  # Get all time series and kymographs
  ratio.data <- GetDataFromRatioKymo(k1, k2, tip.find.par, image.par)
  
  # Analyze all time series
  time.vec <- ratio.data$all.kymo.ts$time
  ts.nms <- c("growth.raw", "growth.smth","tip.fluo", "shank.fluo")
  ts.analysis <- list()
  analysis.tbls <- NULL
  summary.tbls <- NULL
  
  # Save preferences before altering
  do.save.data <- out.par$do.save.data
  out.par$do.save.data <- FALSE
  
  for(var.nm in ts.nms){
    dat <- cbind(time.vec, ratio.data$all.kymo.ts[var.nm])
    colnames(dat) <- c("time", var.nm)
    fl.nm.plt <- paste(toupper(var.nm), "_", fl.nm, sep = "")
    # FUNCTION FROM: 'TimeSeriesAnalysisModule.R'
    if(var.nm %in% c("tip.fluo", "shank.fluo")){
      out.par$ylab <- paste("Fluorescence (", var2.nm, "/", var1.nm, ")", 
                            sep = "")
    }
    ts.analysis[[var.nm]] <- AnalyzeTimeSeries(dat[!is.na(dat[, 2]), ],
                                               ts.par, filter.par, wvlt.par, 
                                               out.par, var.nm, fl.nm.plt)
    # FUNCTION FROM: /src_aux/'AuxFuns.R'
    analysis.tbls <- AddIfNotNull(analysis.tbls, ts.analysis[[var.nm]]$wvlt.tbl)
    summary.tbls <- AddIfNotNull(summary.tbls, ts.analysis[[var.nm]]$summary.stats.tbl)
  }
  if(out.par$do.plot){
    ylm <- c(0, max(ratio.data$all.kymo.ts$tip.loc.smth * image.par$pixel.size) + 5)
    pdf(paste(out.par$out.path, "RatioKymo_", fl.nm, ".pdf", sep = ""), 
        width = 9, height = 6)
    par(mar = c(5, 4.5, 4, 4.5) + 0.1)
    
    # Strongest channel kymograph
    image.plot(x = time.vec, y = 1:dim(k2)[2] * image.par$pixel.size, z = k2, col = tim.colors(256),
          xlab = paste("Time (", image.par$time.unit, ")", sep=""),
          ylab = paste("Length (", image.par$length.unit, ")", sep=""), ylim = ylm,
          main = paste(var2.nm, " kymograph"))
    
    lines(ratio.data$all.kymo.ts$time, ratio.data$all.kymo.ts$tip.loc.smth * image.par$pixel.size, lwd = 1, col ="white", type="o",cex=0.3, pch = 19)
    legend("topleft","Tip location estimate", lwd = 1, col="white",text.col = "white", pch = 19, pt.cex = 0.5,bty="n")#box.col = "white",

    # Weakest channel kymograph
    image.plot(x = time.vec, y = 1:dim(k1)[2] * image.par$pixel.size, z = k1, col = tim.colors(256),
               xlab = paste("Time (", image.par$time.unit, ")", sep=""), ylim = ylm,
               ylab = paste("Length (", image.par$length.unit, ")", sep=""),
               main = paste(var1.nm, " kymograph"))
    
    lines(ratio.data$all.kymo.ts$time, ratio.data$all.kymo.ts$tip.loc.smth * image.par$pixel.size, lwd = 1, col ="white", type="o",cex=0.3, pch = 19)
    legend("topleft","Tip location estimate", lwd = 1, col="white",text.col = "white", pch = 19, pt.cex = 0.5,bty="n")#box.col = "white",
    
    # Tip aligned kymograph
    tip.aligned <- ratio.data$ratio.kymo.lst$tip.aligned
    image.plot(x = time.vec, 
               y = image.par$tip.cut:dim(k2)[2] * image.par$pixel.size,
               z = tip.aligned[, image.par$tip.cut:dim(tip.aligned)[2]], 
               col = tim.colors(256), ylim = ylm,
               xlab = paste("Time (", image.par$time.unit, ")", sep=""),
               ylab = paste("Length (", image.par$length.unit, ")", sep=""),
               main = "Ratiometric tip aligned kymograph")
    ### With region
    image.plot(x = time.vec, 
               y = image.par$tip.cut:dim(k2)[2] * image.par$pixel.size,
               z = tip.aligned[, image.par$tip.cut:dim(tip.aligned)[2]], 
               col = tim.colors(256), ylim = ylm,
               xlab = paste("Time (", image.par$time.unit, ")", sep=""),
               ylab = paste("Length (", image.par$length.unit, ")", sep=""),
               main = "Ratiometric tip aligned kymograph")
    
    rect(xleft = -1, 
         ybottom = (image.par$tip.pixel - (image.par$avg.width / 2)) * image.par$pixel.size,
         xright = max(time.vec)+1,
         ytop = (image.par$tip.pixel + (image.par$avg.width / 2)) * image.par$pixel.size,
         density = 10,
         col = adjustcolor("black", alpha.f = 0.5)
    )
    text(x=max(time.vec)/2,y=(image.par$tip.pixel + (image.par$avg.width / 2)) * image.par$pixel.size + 3,"Tip",cex=1)
    #0+2
    rect(xleft = -1, 
         ybottom = (image.par$shank.pixel - (image.par$avg.width / 2)) * image.par$pixel.size,
         xright = max(time.vec)+1,
         ytop = (image.par$shank.pixel + (image.par$avg.width / 2)) * image.par$pixel.size,
         density = 10,
         col = adjustcolor("black", alpha.f = 0.5)
    )
    text(x=max(time.vec)/2,y=(image.par$shank.pixel + (image.par$avg.width / 2)) * image.par$pixel.size + 3,"Shank",cex=1)
    
    
    dev.off()
  }
  
  if(do.save.data){
    write.csv(analysis.tbls, 
              paste(out.par$out.path, "OscillationAnalysisTbl_", fl.nm, ".csv", 
                    sep = ""), row.names = FALSE)
    write.csv(summary.tbls, 
              paste(out.par$out.path, "SummaryStats_", fl.nm, ".csv", sep = ""), 
              row.names = FALSE)
    
    write.csv(ratio.data$all.kymo.ts, 
              paste(out.par$out.path, "RatioKymoTimeSeries_", fl.nm, ".csv", sep = ""), 
              row.names = FALSE)
    
    write.table(ratio.data$ratio.kymo.lst$tip.aligned, 
                paste(out.par$out.path, "TipAlignedKymo_", fl.nm, ".csv", 
                      sep = ""), row.names = FALSE, col.names = FALSE, 
                sep = ",")
  }
  return(list("ratio.data" = ratio.data,
              "ts.analysis.lst" = ts.analysis,
              "analysis.tbls" = analysis.tbls,
              "summary.tbls" = summary.tbls))
  
}
#-------------------------------------------------------------------------------
GetDataFromRatioKymo <- function(k1, k2, tip.find.par, image.par){
  # Calculate time vector
  time.vec <- 0:(dim(k1)[1] - 1) * image.par$time.step
  
  # Calculate tip location and growth rate
  # FUNCTION FROM: 'TipDetectionModule.R'
  tip.dat <- DetectTipAndGrowthRate(kymo = k2, tip.find.par, image.par)
  
  # Use smooth or raw tip location estimates to generate kymographs
  if(tip.find.par$use.tip.smth){
    tip.loc <- tip.dat$tip.loc.smth
    } else{
    tip.loc <- tip.dat$tip.loc.raw
    }
  
  # Get list of ratiometric kymographs
  # returns "raw", "clean", "tip.aligned" kymographs
  ratio.kymo.lst <- GetRatiometricKymograph(k1 = k1, 
                                            k2 = k2,
                                            tip.loc = tip.loc,
                                            bckgd.qntl = image.par$bckgd.qntl,
                                            tip.mrgn = image.par$tip.mrgn)
  
  # Calculate tip fluorescence time series
  # FUNCTION FROM: /src_kymo/'GetFluoTimeSeries.R'
  tip.fluo.ts <- GetFluoTimeSeries(ratio.kymo.lst$tip.aligned,
                                   ini.ind = image.par$tip.pixel,
                                   avg.width = image.par$avg.width)
  
  # Calculate shank fluorescence time series
  shank.fluo.ts <- GetFluoTimeSeries(ratio.kymo.lst$tip.aligned,
                                     ini.ind = image.par$shank.pixel,
                                     avg.width = image.par$avg.width)
  
  
  
  # Data frame with all time series
  # considering that 'tip.dat' "time", 
  # "tip.loc.raw", "tip.loc.smth", "growth.raw", "growth.smth"
  all.kymo.ts <- cbind.data.frame(tip.dat, # contains time and 4 variables
                                  "tip.fluo" = tip.fluo.ts, 
                                  "shank.fluo" = shank.fluo.ts)
  
  return(list("all.kymo.ts" = all.kymo.ts,
              "ratio.kymo.lst" = ratio.kymo.lst))
  
}
#-------------------------------------------------------------------------------
GetRatiometricKymograph <- function(k1, k2, tip.loc, 
                                    bckgd.qntl = 0.99, 
                                    tip.mrgn = 5){
  
  # Computes the ratiometric kymograph between kymographs of two channels
  #
  # Args:
  #   k1: A matrix to be the numerator, usually weakest channel, e.g. CFP.
  #       must be a kymograph with time increasing from top to bottom with 
  #       background to the right and signal to the left
  #   k2: A matrix to be the denominator, usually strongest channel, e.g. YFP.
  #       must be a kymograph with time increasing from top to bottom with 
  #       background to the right and signal to the left
  #   tip.loc: A vector with tip location estimates (in pixels)
  #   bckgd.qntl: A number corresponding to the quantile of the background 
  #               fluorescence to be removed from each channel
  #   tip.mrgn: The number of pixels from the tip used as margin to calculate 
  #             the background
  # Returns:
  #   A list with ratiometric kymographs containing
  #         raw: Matrix with simply channel 2 divided by channel 1
  #         clean: Matrix of the ratio between channels after background
  #                      subtraction without NaNs and Infs
  #         tip.aligned: Matrix of the clean ratio without any background
  #                            and with the tip aligned to one side
  
  # Subtract background from channel 1 kymograph
  if(is.null(bckgd.qntl)){
      k1.minus.bckgd <- k1
      k2.minus.bckgd <- k2
  } else{
    # FUNCTION FROM: /src_kymo/'SubtractKymoBackground.R'
    k1.minus.bckgd <- SubtractKymoBackground(k1, tip.loc, bckgd.qntl, tip.mrgn)
    # Subtract background from channel 2 kymograph
    k2.minus.bckgd <- SubtractKymoBackground(k2, tip.loc, bckgd.qntl, tip.mrgn)
    }
  
  # Divide kymographs from both channels and clean matrix of NaNs and Infs
  # FUNCTION FROM: /src_kymo/'CleanKymo.R'
  kymo.ratio.clean <- CleanKymo(k2.minus.bckgd / k1.minus.bckgd)
    
  # Align kymograph by the tip and remove background
  # FUNCTION FROM: /src_kymo/'AlignByTip.R'
  kymo.ratio.align <- AlignByTip(tip.loc, kymo.ratio.clean)
  
  kymo.lst.out <- list("raw" = k2 / k1,
                       "clean" = kymo.ratio.clean,
                       "tip.aligned" = kymo.ratio.align)
  
  return(kymo.lst.out)
}
#-------------------------------------------------------------------------------

