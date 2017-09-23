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
#-------------------------------------------------------------------------------
SaveAllTables <- function(pop.tbls, dir.out, pref = ""){
  
  if(!is.null(pop.tbls$summary.stats.tbl)){
    SaveTableInDir(tbl = pop.tbls$summary.stats.tbl,
                   dir = dir.out, 
                   fl.nm = "SummaryStatsTbl",
                   pref = pref)
  }
  
  if(!is.null(pop.tbls$wvlt.tbl)){
    SaveTableInDir(tbl = pop.tbls$wvlt.tbl,
                   dir = dir.out, 
                   fl.nm = "WaveletTable",
                   pref = pref)
  }
  
  
  if(!is.null(pop.tbls$sync.tbl)){
    SaveTableInDir(tbl = pop.tbls$sync,
                   dir = dir.out, 
                   fl.nm = "JointWaveletsTbl",
                   pref = pref)
  }
  return()
}
#-------------------------------------------------------------------------------
SaveTableInDir <- function(tbl, dir, fl.nm, pref = ""){
  write.csv(tbl, paste(dir, pref, fl.nm,".csv", sep = ""))  
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
ReadVpData <- function(fl, data.dir){
  #print(paste(data.dir, fl, sep = ""))
  dat <- read.csv(paste(data.dir, fl, sep = ""), sep=";", stringsAsFactors=FALSE)
  #read.csv2(paste(data.dir, fl, sep = ""), sep=";", stringsAsFactors=FALSE)
  dat <- dat[!is.na(dat[, 2]), ]
  dat[, 1] <- as.numeric(dat[, 1])
  dat[, 2] <- as.numeric(dat[, 2])
  colnames(dat) <- c("time", "flux")
  return(dat)
}
#-------------------------------------------------------------------------------
ReadTrackData <- function(fl, data.dir){
  dat.raw <- read.csv(paste(data.dir, fl, sep = ""))
  dat.raw <- dat.raw[1:length(unique(dat.raw$Image.Plane)), ]
  dat <- cbind(cumsum(dat.raw$Time.Interval) / 60, dat.raw$Velocity * 60)
  if(dat[1,2] == 0){dat <- dat[-1,]}
  dat <- dat[!is.na(dat[, 2]), ]
  dat[, 1] <- as.numeric(dat[, 1])
  dat[, 2] <- as.numeric(dat[, 2])
  colnames(dat) <- c("time", "growth")
  return(dat)
  }
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
ParseFileName <- function(fl.nm){
  nm.vec <- unlist(strsplit(fl.nm, split = "_"))
  nm.vec <- c(nm.vec[1:length(nm.vec) - 1],
              unlist(strsplit(nm.vec[length(nm.vec)], split = "\\.")))
  if(any(nm.vec=="")){nm.vec <- nm.vec[-which(nm.vec=="")]}
  add.fields <- length(nm.vec) - 8
  if(add.fields > 0){
    names(nm.vec) <- c("species", "genotype", "ion", "technique", "date", 
                       "tube", "obs", paste(rep("extra", add.fields, sep = "_"),
                                            1:add.fields), "extension")
  } else{names(nm.vec) <- c(c("species", "genotype", "ion", "technique", "date", 
                              "tube", "obs")[1:(length(nm.vec)-1)], "extension")}
  
  full.id <- paste(nm.vec[-length(nm.vec)], sep = "", collapse = "_")
  id.pref <- paste(nm.vec[1:3], sep = "", collapse = "_")
  id.suff <- paste(nm.vec[5:(length(nm.vec)-1)], sep = "", collapse = "_")
  id.strict.suff <- paste(nm.vec[5:7], sep = "", collapse = "_")
  
  id <- list("full.id" = full.id, "id.pref" = id.pref, "id.suff" = id.suff, 
             "id.strict.suff" = id.strict.suff)
  
  fl.attr <- list("nm.lst" = as.list(nm.vec), "id" = id)
  return(fl.attr)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
JoinTables <- function(pop.tbls, analysis.tbls){
  pop.tbls$summary.stats.tbl <- AddIfNotNull(pop.tbls$summary.stats.tbl, 
                                             analysis.tbls$summary.stats.tbl)
  pop.tbls$wvlt.tbl <- AddIfNotNull(pop.tbls$wvlt.tbl, 
                                    analysis.tbls$wvlt.tbl)
  pop.tbls$sync.tbl <- AddIfNotNull(pop.tbls$sync.tbl, 
                                    analysis.tbls$sync.tbl)
  return(pop.tbls)
}
#-------------------------------------------------------------------------------
