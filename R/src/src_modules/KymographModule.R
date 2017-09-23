#-------------------------------------------------------------------------------
AnalyzeKymograph <- function(kymo, tip.find.par, image.par, ts.par, filter.par, 
                             wvlt.par, out.par, fluo.var.nm, fl.nm = ""){
  # Get all time series and kymographs
  kymo.data <- GetDataFromKymo(kymo, tip.find.par, image.par, filter.par)
  # Analyze all time series
  time.vec <- kymo.data$all.kymo.ts$time
  ts.nms <- c("growth.raw", "growth.smth","tip.fluo", "shank.fluo")
  ts.analysis <- list()
  analysis.tbls <- NULL
  summary.tbls <- NULL
  
  # Save preferences before altering
  do.save.data <- out.par$do.save.data
  out.par$do.save.data <- FALSE
  
  for(var.nm in ts.nms){
    dat <- cbind(time.vec, kymo.data$all.kymo.ts[var.nm])
    colnames(dat) <- c("time", var.nm)
    fl.nm.plt <- paste(toupper(var.nm), "_", fl.nm, sep = "")
    # FUNCTION FROM: 'TimeSeriesAnalysisModule.R'
    if(var.nm %in% c("tip.fluo", "shank.fluo")){
      out.par$ylab <- paste(fluo.var.nm, " fluorescence (AU)", sep = "")
    }
    ts.analysis[[var.nm]] <- AnalyzeTimeSeries(dat[!is.na(dat[, 2]), ],
                                               ts.par, filter.par, wvlt.par, 
                                               out.par, var.nm, fl.nm.plt)
    # FUNCTION FROM: /src_aux/'AuxFuns.R'
    analysis.tbls <- AddIfNotNull(analysis.tbls, ts.analysis[[var.nm]]$wvlt.tbl)
    summary.tbls <- AddIfNotNull(summary.tbls, ts.analysis[[var.nm]]$summary.stats.tbl)
  }
  
  if(out.par$do.plot){
    ylm <- c(0, max(kymo.data$all.kymo.ts$tip.loc.smth * image.par$pixel.size) + 5)
    pdf(paste(out.par$out.path, "Kymo_", fl.nm, ".pdf", sep = ""), 
        width = 9, height = 6)
    par(mar = c(5, 4.5, 4, 4.5) + 0.1)
    
    # Original kymograph
    image.plot(x = time.vec, y = 1:dim(kymo)[2] * image.par$pixel.size, z = kymo, col = tim.colors(256),
               ylim = ylm,
               xlab = paste("Time (", image.par$time.unit, ")",sep = ""),
               ylab = paste("Length (", image.par$length.unit, ")",sep = ""),
               main = paste(fluo.var.nm, " kymograph"))
    
    lines(kymo.data$all.kymo.ts$time, kymo.data$all.kymo.ts$tip.loc.smth * image.par$pixel.size, lwd = 1, col ="white", type="o",cex=0.3, pch = 19)
    legend("topleft","Tip location estimate", lwd = 1, col="white",text.col = "white", pch = 19, pt.cex = 0.5,bty="n")#box.col = "white",
    
    # Tip aligned kymograph
    tip.aligned <- kymo.data$kymo.lst$tip.aligned
    cut.time <- c(1:image.par$time.cut, 
                  (dim(tip.aligned)[1] - image.par$time.cut):dim(tip.aligned)[1])
    image.plot(x = time.vec[-cut.time], 
               y = 1:dim(tip.aligned)[2] * image.par$pixel.size,
               z = tip.aligned[-cut.time, ], 
               col = tim.colors(256), ylim = ylm,
               xlab = paste("Time (", image.par$time.unit, ")",sep = ""),
               ylab = paste("Length (", image.par$length.unit, ")",sep = ""),
               main = "Filtered tip aligned kymograph")
    
    # Tip aligned kymograph WITH ROIs
    image.plot(x = time.vec[-cut.time], 
               y = 1:dim(tip.aligned)[2] * image.par$pixel.size,
               z = tip.aligned[-cut.time, ], 
               col = tim.colors(256), ylim = ylm,
               xlab = paste("Time (", image.par$time.unit, ")",sep = ""),
               ylab = paste("Length (", image.par$length.unit, ")",sep = ""),
               main = "Filtered tip aligned kymograph")
    
    rect(xleft = -1, 
         ybottom = (image.par$tip.pixel - (image.par$avg.width / 2)) * image.par$pixel.size,
         xright = max(time.vec)+1,
         ytop = (image.par$tip.pixel + (image.par$avg.width / 2)) * image.par$pixel.size,
         density = 10,
         col = adjustcolor("black", alpha.f = 0.5)
    )
    text(x= max(time.vec)/2, y=(image.par$tip.pixel + (image.par$avg.width / 2)) * image.par$pixel.size + 3,"Tip",cex=1)
    #0+0.5
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
    
    write.csv(kymo.data$all.kymo.ts, 
              paste(out.par$out.path, "KymoTimeSeries_", fl.nm, ".csv", sep = ""), 
              row.names = FALSE)
    
    write.table(kymo.data$kymo.lst$tip.aligned, 
                paste(out.par$out.path, "TipAlignedKymo_", fl.nm, ".csv", 
                      sep = ""), row.names = FALSE, col.names = FALSE, 
                sep = ",")
  }
  
  return(list("kymo.data" = kymo.data,
              "ts.analysis.lst" = ts.analysis,
              "wvlt.tbl" = analysis.tbls,
              "summary.stats.tbl" = summary.tbls))
  
}
#-------------------------------------------------------------------------------
GetDataFromKymo <- function(kymo, tip.find.par, image.par, filter.par){
  # Calculate time vector
  time.vec <- 0:(dim(kymo)[1] - 1) * image.par$time.step
  time.step <- GetSamplingIntervalMode(time.vec, dg = 6)
  # Calculate sampling frequency
  sampling.freq <- 1 / time.step
  if(is.null(filter.par$smth.vec)){
    filter.par$smth.vec <- GetSmoothVec(sampling.freq,
                                    series.length = length(time.vec) * time.step,
                                    low.per = filter.par$low.per, 
                                    high.per = filter.par$high.per)
  } else if(filter.par$smth.vec == "skip"){ 
    filter.par$smth.vec <- GetSmoothVec(sampling.freq,
                                        series.length = length(time.vec) * time.step,
                                        low.per = filter.par$low.per,
                                        high.per = filter.par$high.per)}
  
  # Calculate tip location and growth rate
  # FUNCTION FROM: 'TipDetectionModule.R'
  tip.dat <- DetectTipAndGrowthRate(kymo, tip.find.par, image.par)
  
  # Use smooth or raw tip location estimates to generate kymographs
  if(tip.find.par$use.tip.smth){
    tip.loc <- tip.dat$tip.loc.smth
  } else{
    tip.loc <- tip.dat$tip.loc.raw
  }
  
  kymo.lst <- GetKymograph(kymo = kymo, tip.loc = tip.loc, 
                           filter.par$filter.type, filter.par$smth.vec)
  
  # Calculate tip fluorescence time series
  # FUNCTION FROM: /src_kymo/'GetFluoTimeSeries.R'
  tip.fluo.ts <- GetFluoTimeSeries(kymo.lst$tip.aligned,
                                   ini.ind = image.par$tip.pixel,
                                   avg.width = image.par$avg.width)
  
  # Calculate shank fluorescence time series
  shank.fluo.ts <- GetFluoTimeSeries(kymo.lst$tip.aligned,
                                     ini.ind = image.par$shank.pixel,
                                     avg.width = image.par$avg.width)
  
  
  
  # Data frame with all time series
  # considering that 'tip.dat' "time", 
  # "tip.loc.raw", "tip.loc.smth", "growth.raw", "growth.smth"
  all.kymo.ts <- cbind.data.frame(tip.dat, # contains time and 4 variables
                                  "tip.fluo" = tip.fluo.ts, 
                                  "shank.fluo" = shank.fluo.ts)
  
  return(list("all.kymo.ts" = all.kymo.ts,
              "kymo.lst" = kymo.lst))
  
}
#-------------------------------------------------------------------------------
# TODO: USE MODULE TO GET SERIES THROUGHOUT THE TUBE
# TODO: COMPARE WITH POOR MAN'S RATIOMETRIC DeltaF/F0
GetKymograph <- function(kymo, tip.loc, filter.type = NULL, smth.vec = c(F,T,T,T,F,F)){
  
  if(is.null(filter.type)){
    kymo.filt <- kymo
  } else if(filter.type == "Filter2DWavelet"){
    kymo.filt <- Filter2DWavelet(imaj = kymo, smth.vec = smth.vec)
  }
  kymo.align <- AlignByTip(tip.loc, kymo.filt)
  
  kymo.lst.out <- list("raw" = kymo,
                       "filt" = kymo.filt,
                       "tip.aligned" = kymo.align)
  
  return(kymo.lst.out)
}
#-------------------------------------------------------------------------------