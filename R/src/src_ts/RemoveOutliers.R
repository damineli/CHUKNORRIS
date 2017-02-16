#-------------------------------------------------------------------------------
RemoveOutliers <- function(dat, sampling.freq, outlier.removal = "tso",
                           integrate.ts = TRUE, out.thrsh = NULL, smth.n = 5){
  if(outlier.removal == "tso"){
    y.old <- dat[, 2]
    if(integrate.ts){y.old <- cumsum(dat[, 2])}
    tso.mdl <- tsoutliers::tso(ts(y.old, frequency =  sampling.freq), 
                               types = c("LS", "TC")) # to try harder add:
                                # ",maxit = 100, maxit.iloop = 10" to try harder
    y.new <- tso.mdl$yadj
    if(integrate.ts){y.new <- c(y.new[1], diff(y.new))}
  } else if(outlier.removal == "mad"){
    mod.z <- (0.6745 * (diff(dat[,2 ]) - median(diff(dat[, 2])))) / 
      mad(diff(dat[, 2])) #(modified Z score Iglewicz and Hoaglin, 1993)
    outlier.indx <- which(abs(mod.z) > 3.5)
    y.new <- FillOutliers(dat, outlier.indx, smth.n = smth.n)[, 2]
  } else if(outlier.removal == "thrsh"){
    outlier.indx <- which(abs(dat[, 2]) > out.thrsh)
    y.new <- FillOutliers(dat, outlier.indx, smth.n = smth.n)[, 2]
  } else if(outlier.removal == "iqr"){
    qntls <- quantile(dat[, 2], probs = c(0.25, 0.75))
    upper.out <- which(dat[, 2] > (qntls[2] + 3 * diff(qntls)))
    lower.out <- which(dat[, 2] < (qntls[1] - 3 * diff(qntls)))
    outlier.indx <- c(upper.out, lower.out)
    y.new <- FillOutliers(dat, outlier.indx, smth.n = smth.n)[, 2]
  } else{y.new <- dat[, 2]}
  dat.clean <- cbind(dat[, 1], y.new)
  colnames(dat.clean) <- colnames(dat)
  return(dat.clean)
}

#-------------------------------------------------------------------------------
FillOutliers <- function(dat, outlier.indx, smth.n = 5){
  if(length(outlier.indx) > 0){
    dat.add <- InterpolateWithLoess(dat[-outlier.indx, ], smth.n = smth.n,
                                    interpol.times = dat[outlier.indx, 1])
    dat.interp <- dat
    dat.interp[outlier.indx, 2] <- dat.add[, 2]
  } else{dat.interp <- dat}
  return(dat.interp)
}
#-------------------------------------------------------------------------------
